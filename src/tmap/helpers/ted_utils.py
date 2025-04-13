"""
Utility functions for Tree Edit Distance (TED) calculations using the apted package.
This module provides functions to convert between tmap's SynthesisTree format and the
format required by the apted package, as well as utility functions for computing
pairwise distances between trees.
"""

from typing import List, Dict, Any, Union, Tuple, Optional
import numpy as np
from tmap.core import SynthesisTree, SynthesisTreeNode

try:
    from apted import APTED, Config
    APTED_AVAILABLE = True
except ImportError:
    APTED_AVAILABLE = False

def get_reaction_classifications(classification: str) -> list[int]:
    """Splits AiZynthFinder matadata into numbers"""

    rc = classification.split(" ")[0]  # Take numeric classes, ignore text name.
    rcn = [int(i) for i in rc.split(".")]
    return rcn


class CustomConfig(Config):
    """Custom APTED config for CAZP"""

    def delete(self, node):
        """Calculates the cost of deleting a node"""
        return 1

    def insert(self, node):
        """Calculates the cost of inserting a node"""
        return 1

    def rename(self, node1, node2):
        if node1["type"] != node2["type"]:
            return 1

        if node1["type"] == "reaction":

            rc1n = get_reaction_classifications(node1["metadata"]["classification"])
            rc2n = get_reaction_classifications(node2["metadata"]["classification"])
            if len(rc1n) != len(rc2n):
                return 1
            diffs = ~np.equal(rc1n, rc2n)
            if diffs[0]:  # Diff in first number.
                dist = 0.8
            elif len(diffs) > 1 and diffs[1]:
                dist = 0.5
            elif len(diffs) > 2 and diffs[2]:
                dist = 0
            else:  # No diff.
                dist = 0
            return dist

        if node1["type"] == "mol":
            return 0

    def children(self, node):
        """Get children of a tree"""
        return node.get("children", [])


def compute_ted_with_apted(tree1: Dict, tree2: Dict) -> float:
    """
    Compute Tree Edit Distance using the APTED algorithm with custom config.
    
    Args:
        tree1: First tree in dict format
        tree2: Second tree in dict format
        
    Returns:
        float: The Tree Edit Distance between the trees
    """
    if not APTED_AVAILABLE:
        raise ImportError("APTED package is required for Tree Edit Distance calculation. Please install it with 'pip install apted'.")
    
    apted = APTED(tree1, tree2, CustomConfig())
    return apted.compute_edit_distance()

def convert_to_apted_tree(tree: Union[SynthesisTree, SynthesisTreeNode]) -> Dict:
    """
    Convert a SynthesisTree or SynthesisTreeNode to the format required by the apted package.
    
    Args:
        tree: The SynthesisTree or SynthesisTreeNode to convert
        
    Returns:
        Dict: The tree in apted format
    """
    if isinstance(tree, SynthesisTree):
        # If we received a SynthesisTree, use its root
        node = tree.root
    else:
        # Otherwise, assume it's a SynthesisTreeNode
        node = tree
    
    # Determine node type from data or default to "mol"
    node_type = node.data.get("type", "mol")
    
    # Create a node dictionary compatible with the custom APTED config
    apted_node = {
        "type": node_type,
        "label": node.label,
        "metadata": node.data.get("metadata", {}),
        "children": []
    }
    
    # Recursively convert children
    if node.children:
        apted_node["children"] = [convert_to_apted_tree(child) for child in node.children]
    
    return apted_node

def compute_ted(tree1: Union[SynthesisTree, SynthesisTreeNode], 
               tree2: Union[SynthesisTree, SynthesisTreeNode]) -> float:
    """
    Compute the Tree Edit Distance between two trees.
    
    Args:
        tree1: First tree
        tree2: Second tree
        
    Returns:
        float: The Tree Edit Distance between the trees
    """
    if APTED_AVAILABLE:
        # Convert trees to APTED format
        apted_tree1 = convert_to_apted_tree(tree1)
        apted_tree2 = convert_to_apted_tree(tree2)
        
        # Compute TED using APTED
        return compute_ted_with_apted(apted_tree1, apted_tree2)
    else:
        # Fallback to a simple approximation if APTED is not available
        def count_nodes(node: SynthesisTreeNode) -> int:
            return 1 + sum(count_nodes(child) for child in node.children)
        
        if isinstance(tree1, SynthesisTree):
            node1 = tree1.root
        else:
            node1 = tree1
            
        if isinstance(tree2, SynthesisTree):
            node2 = tree2.root
        else:
            node2 = tree2
        
        count1 = count_nodes(node1)
        count2 = count_nodes(node2)
        
        # Crude approximation based on node count difference
        return abs(count1 - count2)

def compute_pairwise_distances(trees: List[Union[SynthesisTree, SynthesisTreeNode]]) -> np.ndarray:
    """
    Compute pairwise TED distances between all trees in the given list.
    
    Args:
        trees: List of trees
        
    Returns:
        np.ndarray: n x n matrix of pairwise distances
    """
    n = len(trees)
    distances = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            dist = compute_ted(trees[i], trees[j])
            distances[i, j] = dist
            distances[j, i] = dist
    
    return distances

def get_k_nearest_neighbors(distances: np.ndarray, k: int) -> List[Tuple[int, int, float]]:
    """
    Get the k nearest neighbors for each node based on pairwise distances.
    
    Args:
        distances: n x n matrix of pairwise distances
        k: Number of nearest neighbors to find for each node
        
    Returns:
        List[Tuple[int, int, float]]: List of (source, target, distance) tuples
    """
    n = distances.shape[0]
    edge_list = []
    
    for i in range(n):
        # Get indices of k nearest neighbors for node i
        # We exclude the distance to self (which is 0)
        neighbor_indices = np.argsort(distances[i])
        # Skip the first one (which is the node itself with distance 0)
        neighbor_indices = neighbor_indices[1:k+1]
        
        for j in neighbor_indices:
            edge_list.append((i, j, distances[i, j]))
    
    return edge_list 