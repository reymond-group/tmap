"""
Example of using the TEDLayoutGenerator to visualize synthesis route trees.
This example creates a few sample synthesis trees and visualizes them using
the Tree Edit Distance (TED) to group similar trees together.
"""
import examples.synthesis_routes.use_local_tmap as use_local_tmap
import tmap as tm
import numpy as np
from matplotlib import pyplot as plt
from tmap.core import SynthesisTree, SynthesisTreeNode
from tmap.layout_generators import TEDLayoutGenerator
import random

def create_sample_trees(n_trees: int = 10, max_depth: int = 3, max_children: int = 3) -> list:
    """
    Create a list of sample synthesis trees for demonstration purposes.
    
    Args:
        n_trees: Number of trees to create
        max_depth: Maximum depth of each tree
        max_children: Maximum number of children for each node
        
    Returns:
        list: List of SynthesisTree objects
    """
    trees = []
    
    for i in range(n_trees):
        # Create a root node for the tree
        root = SynthesisTreeNode(label=f"Product {i}")
        
        # Function to recursively create child nodes
        def add_children(node, depth):
            if depth >= max_depth:
                return
                
            # Randomly determine the number of children
            n_children = random.randint(0, max_children)
            
            for j in range(n_children):
                child = SynthesisTreeNode(label=f"Intermediate {i}.{depth}.{j}")
                node.add_child(child)
                
                # Recursively add children to this child
                add_children(child, depth + 1)
        
        # Add children to the root node
        add_children(root, 1)
        
        # Create a tree from the root node
        tree = SynthesisTree(root)
        trees.append(tree)
    
    return trees

def main():
    """
    Main function to demonstrate the TEDLayoutGenerator.
    """
    # Create a list of sample synthesis trees
    n_trees = 20
    print(f"Creating {n_trees} sample synthesis trees...")
    trees = create_sample_trees(n_trees, max_depth=3, max_children=2)
    
    # Print some information about the trees
    for i, tree in enumerate(trees):
        def count_nodes(node):
            return 1 + sum(count_nodes(child) for child in node.children)
        
        n_nodes = count_nodes(tree.root)
        print(f"Tree {i}: {n_nodes} nodes")
    
    # Create a TED layout generator
    print("Creating TED layout generator...")
    layout_generator = TEDLayoutGenerator(
        k=5,  # Number of nearest neighbors
        create_mst=True,  # Create minimum spanning tree
        keep_knn=True,    # Keep k-nearest neighbors graph
    )
    
    # Generate the layout
    print("Generating layout...")
    embedding = layout_generator.layout(trees)
    
    # Extract coordinates and edges
    x, y = embedding.get_coordinates()
    s, t, w = embedding.get_edges()
    
    # Plot the layout
    print("Plotting layout...")
    plt.figure(figsize=(10, 10))
    
    # Plot edges with width proportional to weight (inverse of distance)
    for i in range(len(s)):
        plt.plot(
            [x[s[i]], x[t[i]]],
            [y[s[i]], y[t[i]]],
            "k-",
            linewidth=0.5,
            alpha=0.5,
            zorder=1,
        )
    
    # Plot vertices with size proportional to number of nodes in the tree
    sizes = []
    for tree in trees:
        def count_nodes(node):
            return 1 + sum(count_nodes(child) for child in node.children)
        
        sizes.append(50 * count_nodes(tree.root))
    
    plt.scatter(x, y, s=sizes, c=range(len(trees)), cmap="viridis", zorder=2)
    
    # Add labels for each point
    for i in range(len(trees)):
        plt.annotate(
            f"{i}",
            (x[i], y[i]),
            textcoords="offset points",
            xytext=(0, 5),
            ha="center",
        )
    
    plt.title("Synthesis Route Trees Grouped by Tree Edit Distance")
    plt.tight_layout()
    
    # Save the figure
    plt.savefig("synthesis_trees_ted.png")
    print("Figure saved as 'synthesis_trees_ted.png'")
    
    # Show the figure
    plt.show()

if __name__ == "__main__":
    main() 