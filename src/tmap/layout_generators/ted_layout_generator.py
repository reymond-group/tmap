from typing import Iterable, Callable, List, Any, Dict, Tuple, Union
import tmap as tm
import numpy as np
from tmap.core import TMAPEmbedding, SynthesisTree, SynthesisTreeNode
from tmap.helpers.ted_utils import compute_ted, compute_pairwise_distances, get_k_nearest_neighbors
from .base_layout_generator import BaseLayoutGenerator

class TEDLayoutGenerator(BaseLayoutGenerator):
    """
    Layout generator using Tree Edit Distance (TED) for measuring similarities between
    synthesis route trees. This class uses the apted package to compute the
    pairwise distances between tree structures.
    """
    def __init__(
        self,
        distance_function: Callable[[Any, Any], float] = None,
        create_mst: bool = True,
        keep_knn: bool = False,
        k: int = 10,
        kc: int = 10,
        fme_iterations: int = 100,
        fme_threads: int = 4,
        fme_precision: int = 4,
        sl_repeats: int = 1,
        sl_extra_scaling_steps: int = 2,
        sl_scaling_min: float = 1,
        sl_scaling_max: float = 1,
        sl_scaling_type: tm.ScalingType = tm.ScalingType.RelativeToDrawing,
        mmm_repeats: int = 1,
        placer: tm.Placer = tm.Placer.Barycenter,
        merger: tm.Merger = tm.LocalBiconnected,
        merger_factor: float = 2,
        merger_adjustment: int = 0,
        node_size: float = 1 / 65,
        ted_threshold: float = 5.0,
        force_connect: bool = False,
    ) -> None:
        """
        Initialize the TED layout generator.
        
        Args:
            distance_function: Function to compute TED between two trees. If None, a default function
                              will be implemented using the apted package (to be implemented later).
            create_mst: Whether to create a minimum spanning tree from the k-nearest neighbors graph.
            keep_knn: Whether to keep the k-nearest neighbors graph edges in addition to the MST.
            k: Number of nearest neighbors to use.
            kc: Number of candidates to consider for each query point.
            fme_iterations: Number of iterations for the Fast Multipole Embedder.
            fme_threads: Number of threads to use for the Fast Multipole Embedder.
            fme_precision: Precision for the Fast Multipole Embedder.
            sl_repeats: Number of times to repeat the scaling layout algorithm.
            sl_extra_scaling_steps: Number of additional scaling steps.
            sl_scaling_min: Minimum scaling factor.
            sl_scaling_max: Maximum scaling factor.
            sl_scaling_type: Type of scaling to use.
            mmm_repeats: Number of times to repeat the multilevel layout algorithm.
            placer: Strategy for placing vertices in the layout.
            merger: Strategy for merging vertices during layout.
            merger_factor: Parameter for the merger strategy.
            merger_adjustment: Adjustment for the merger.
            node_size: Size of nodes in the layout.
            ted_threshold: Maximum Tree Edit Distance (TED) allowed for connecting nodes.
                          Trees with TED less than or equal to this threshold will be connected.
                          Lower values = more strict similarity requirement. Default: 5.0
            force_connect: Whether to force a connected graph even if it breaks the threshold.
                          If True, will create connections to ensure all nodes are in one graph.
                          Default: False.
        """
        super().__init__(
            create_mst,
            keep_knn,
            k,
            kc,
            fme_iterations,
            fme_threads,
            fme_precision,
            sl_repeats,
            sl_extra_scaling_steps,
            sl_scaling_min,
            sl_scaling_max,
            sl_scaling_type,
            mmm_repeats,
            placer,
            merger,
            merger_factor,
            merger_adjustment,
            node_size,
        )

        # Store the distance function or use default when implemented
        self.distance_function = distance_function
        self.ted_threshold = ted_threshold
        self.force_connect = force_connect
        
    def _compute_ted_distance(self, tree1: Union[SynthesisTree, SynthesisTreeNode], 
                             tree2: Union[SynthesisTree, SynthesisTreeNode]) -> float:
        """
        Compute the Tree Edit Distance between two trees.
        
        Args:
            tree1: First tree
            tree2: Second tree
            
        Returns:
            float: The Tree Edit Distance between the trees
        """
        if self.distance_function is not None:
            return self.distance_function(tree1, tree2)
        
        # Use the helper function from ted_utils
        return compute_ted(tree1, tree2)

    def _find_connected_components(self, edge_list, n):
        """
        Find connected components in the graph defined by edge_list.
        
        Args:
            edge_list: List of edges in the format (source, target, weight)
            n: Number of nodes
            
        Returns:
            List of connected components, where each component is a list of node indices
        """
        # Create adjacency list
        adj_list = [[] for _ in range(n)]
        for s, t, _ in edge_list:
            adj_list[s].append(t)
            adj_list[t].append(s)
        
        # Find connected components using DFS
        visited = [False] * n
        components = []
        
        def dfs(node, component):
            visited[node] = True
            component.append(node)
            for neighbor in adj_list[node]:
                if not visited[neighbor]:
                    dfs(neighbor, component)
        
        for i in range(n):
            if not visited[i]:
                component = []
                dfs(i, component)
                components.append(component)
        
        return components
    
    def layout_from_edge_list(self, n: int, edge_list: List[Tuple[int, int, float]], 
                               create_mst: bool = None, keep_knn: bool = None) -> TMAPEmbedding:
        """
        Override the base class method to prevent it from grouping components that don't meet
        the similarity threshold.
        
        Args:
            n: Number of nodes
            edge_list: List of edges as (source, target, weight) tuples
            create_mst: Whether to create a minimum spanning tree
            keep_knn: Whether to keep the k-nearest neighbors graph
            
        Returns:
            TMAPEmbedding: The embedding
        """
        # Just pass through to parent implementation
        # We already processed the edge list to respect TED thresholds
        return super().layout_from_edge_list(n, edge_list, create_mst, keep_knn)

    def layout(self, trees: Iterable[Union[SynthesisTree, SynthesisTreeNode]], 
              create_mst: bool = None, keep_knn: bool = None) -> TMAPEmbedding:
        """
        Generate a layout for the given trees using Tree Edit Distance (TED).
        
        Args:
            trees: Iterable of trees to lay out
            create_mst: Whether to create a minimum spanning tree
            keep_knn: Whether to keep the k-nearest neighbors graph edges
            
        Returns:
            TMAPEmbedding: The generated layout
        """
        if create_mst is None:
            create_mst = self.create_mst
        if keep_knn is None:
            keep_knn = self.keep_knn
            
        # Convert trees to list for random access
        tree_list = list(trees)
        n = len(tree_list)
        
        if n == 0:
            raise ValueError("No trees provided for layout")
            
        # Compute pairwise distances using TED
        distances = compute_pairwise_distances(tree_list)
        
        # Print raw TED distance statistics for diagnostics
        print(f"TED distance statistics:")
        print(f"  min={np.min(distances):.4f}, max={np.max(distances):.4f}, mean={np.mean(distances):.4f}, median={np.median(distances):.4f}")
        print(f"  percentiles: 10%={np.percentile(distances, 10):.4f}, 25%={np.percentile(distances, 25):.4f}, 50%={np.percentile(distances, 50):.4f}, 75%={np.percentile(distances, 75):.4f}, 90%={np.percentile(distances, 90):.4f}")
        print(f"Applying TED threshold: {self.ted_threshold}")
        
        # Check how many pairs would meet the threshold
        ted_below_threshold = np.sum(distances <= self.ted_threshold) - n  # Subtract diagonal (self-comparisons)
        total_possible_edges = (n * (n-1)) // 2
        print(f"Found {ted_below_threshold} TED values below threshold {self.ted_threshold} out of {total_possible_edges} possible pairs")
        print(f"This represents {ted_below_threshold / total_possible_edges * 100:.2f}% of possible edges")
        
        # Create edge list strictly applying the TED threshold
        edge_list = []
        edges_added = 0
        
        # Only connect nodes if they meet the threshold
        for i in range(n):
            for j in range(i+1, n):  # Avoid duplicate edges
                if distances[i, j] <= self.ted_threshold:
                    edge_list.append((i, j, distances[i, j]))
                    edges_added += 1
        
        print(f"Created {edges_added} edges that meet the TED threshold {self.ted_threshold}")
        
        # Find connected components
        components = self._find_connected_components(edge_list, n)
        num_components = len(components)
        print(f"Number of connected components: {num_components}")
        if num_components > 1:
            for i, comp in enumerate(components):
                print(f"  Component {i+1}: {len(comp)} nodes")
        else:
            print(f"  Single component with {len(components[0])} nodes")
        
        # If force_connect is enabled and we have multiple components, add minimum connections
        if self.force_connect and num_components > 1:
            print("Forcing connection between components (may break TED threshold)")
            # Find the best connections between components
            for i in range(num_components):
                for j in range(i+1, num_components):
                    # Find the closest pair of nodes between components i and j
                    best_ted = float('inf')
                    best_edge = None
                    for ni in components[i]:
                        for nj in components[j]:
                            ted = distances[ni, nj]
                            if ted < best_ted:
                                best_ted = ted
                                best_edge = (ni, nj, ted)
                    
                    if best_edge:
                        edge_list.append(best_edge)
                        print(f"  Added connection between components {i+1} and {j+1} with TED {best_ted:.4f}")
            
            # Verify we now have a single connected component
            components = self._find_connected_components(edge_list, n)
            print(f"After forcing connection: {len(components)} connected components")
        
        # IMPORTANT: We need to disable MST creation in the base class if we have multiple components
        # and force_connect is False, otherwise it will connect dissimilar components
        base_create_mst = create_mst
        if num_components > 1 and not self.force_connect and create_mst:
            print(f"WARNING: Disabling MST creation in base class to preserve {num_components} disconnected components")
            base_create_mst = False
            
            # We'll create MST within each component separately ourselves
            mst_edge_list = []
            for component in components:
                if len(component) <= 1:
                    continue  # Skip singleton components
                    
                # Create subgraph for this component
                subgraph_edges = [(edge[0], edge[1], edge[2]) for edge in edge_list 
                                 if edge[0] in component and edge[1] in component]
                
                # Create MST for this component using a simple implementation
                # This assumes the component is connected, which it should be
                component_nodes = set(component)
                mst_nodes = {component[0]}  # Start with first node
                while len(mst_nodes) < len(component_nodes):
                    best_edge = None
                    best_weight = float('inf')
                    
                    # Find the best edge connecting a node in the MST to a node outside
                    for u, v, w in subgraph_edges:
                        if (u in mst_nodes and v not in mst_nodes) or (v in mst_nodes and u not in mst_nodes):
                            if w < best_weight:
                                best_edge = (u, v, w)
                                best_weight = w
                    
                    if best_edge:
                        mst_edge_list.append(best_edge)
                        mst_nodes.add(best_edge[0])
                        mst_nodes.add(best_edge[1])
                    else:
                        # If we can't find an edge, the component might not be fully connected
                        # Just use all edges for this component
                        mst_edge_list.extend(subgraph_edges)
                        break
            
            # Use our component-wise MSTs
            edge_list = mst_edge_list
            print(f"Created {len(edge_list)} MST edges across {num_components} components")
        
        # Use the base class to create the layout from the edge list
        return super().layout_from_edge_list(n, edge_list, base_create_mst, keep_knn)