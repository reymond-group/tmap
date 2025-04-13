"""
Example of using the TEDLayoutGenerator to visualize synthesis route trees from AiZynthFinder output.
This script loads trees from an AiZynthFinder output.json.gz file and visualizes them
using Tree Edit Distance (TED) to group similar synthesis routes together.
"""
import use_local_tmap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from tmap.core.synthesis_tree import SynthesisTree, SynthesisTreeNode
from tmap.layout_generators.ted_layout_generator import TEDLayoutGenerator
from tmap.helpers.ted_utils import compute_ted
import argparse
import os
import json
import gzip

def convert_aizynthfinder_node_to_synthesis_node(node_data, parent_label=""):
    """
    Convert an AiZynthFinder node to a SynthesisTreeNode.
    
    Args:
        node_data: Node data from AiZynthFinder
        parent_label: Label of the parent node (for generating unique labels)
        
    Returns:
        SynthesisTreeNode: Converted node
    """
    # Extract node metadata
    node_type = node_data.get("type", "mol")
    
    # Create a label based on the node type and data
    if node_type == "mol":
        # For molecules, use SMILES if available, otherwise a generic label
        label = node_data.get("smiles", f"Molecule-{parent_label}")
    elif node_type == "reaction":
        # For reactions, use classification if available
        classification = node_data.get("metadata", {}).get("classification", "")
        label = f"Reaction-{classification}-{parent_label}"
    else:
        # Generic label for other node types
        label = f"{node_type}-{parent_label}"
    
    # Create a SynthesisTreeNode with the extracted data
    node = SynthesisTreeNode(
        label=label,
        data={
            "type": node_type,
            "metadata": node_data.get("metadata", {})
        }
    )
    
    # Recursively convert children
    if "children" in node_data and node_data["children"]:
        for i, child_data in enumerate(node_data["children"]):
            child_node = convert_aizynthfinder_node_to_synthesis_node(
                child_data, 
                parent_label=f"{parent_label}-{i}"
            )
            node.add_child(child_node)
    
    return node

def convert_aizynthfinder_tree_to_synthesis_tree(tree_data, tree_idx=0):
    """
    Convert an AiZynthFinder tree to a SynthesisTree.
    
    Args:
        tree_data: Tree data from AiZynthFinder
        tree_idx: Index of the tree (for generating unique labels)
        
    Returns:
        SynthesisTree: Converted tree
    """
    # Convert the root node
    root_node = convert_aizynthfinder_node_to_synthesis_node(tree_data, parent_label=f"tree{tree_idx}")
    
    # Create a SynthesisTree with the converted root
    return SynthesisTree(root_node)

def load_aizynthfinder_trees(filepath):
    """
    Load trees from an AiZynthFinder output file.
    
    Args:
        filepath: Path to the AiZynthFinder output.json.gz file
        
    Returns:
        list: List of SynthesisTree objects
    """
    # Load the JSON data from the gzipped file
    try:
        data = pd.read_json(filepath, orient="table")
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        # Try alternate loading method
        try:
            with gzip.open(filepath, 'rt') as f:
                data = json.load(f)
                if isinstance(data, dict) and "trees" in data:
                    data = pd.DataFrame({"trees": data["trees"]})
                else:
                    data = pd.DataFrame({"trees": [data]})
        except Exception as e:
            print(f"Error with alternate loading method: {e}")
            return []
    
    synthesis_trees = []
    
    # Extract trees and convert them to SynthesisTree objects
    if hasattr(data, "trees"):
        # Iterate through all molecules in the batch
        for mol_idx, mol_trees in enumerate(data.trees):
            # Iterate through all trees for this molecule
            for tree_idx, tree in enumerate(mol_trees):
                synthesis_tree = convert_aizynthfinder_tree_to_synthesis_tree(tree, tree_idx=f"{mol_idx}_{tree_idx}")
                synthesis_trees.append(synthesis_tree)
                
    return synthesis_trees

def main():
    """
    Main function to demonstrate the TEDLayoutGenerator with AiZynthFinder trees.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Visualize AiZynthFinder trees using Tree Edit Distance")
    parser.add_argument("input_file", help="Path to the AiZynthFinder output.json.gz file")
    parser.add_argument("--k", type=int, default=5, help="Number of nearest neighbors (default: 5)")
    parser.add_argument("--output", default="aizynthfinder_trees_ted.png", help="Output image file (default: aizynthfinder_trees_ted.png)")
    args = parser.parse_args()
    
    # Check if the input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        return
    
    # Load the trees from the AiZynthFinder output file
    print(f"Loading trees from {args.input_file}...")
    trees = load_aizynthfinder_trees(args.input_file)
    
    if not trees:
        print("No trees found. Exiting.")
        return
    
    print(f"Loaded {len(trees)} trees.")
    
    # Create a TED layout generator
    print("Creating TED layout generator...")
    layout_generator = TEDLayoutGenerator(
        k=args.k,       # Number of nearest neighbors
        create_mst=True, # Create minimum spanning tree
        keep_knn=True,   # Keep k-nearest neighbors graph
    )
    
    # Generate the layout
    print("Generating layout...")
    embedding = layout_generator.layout(trees)
    
    # Extract coordinates and edges from embedding properties
    # The TMAPEmbedding object provides x and y as properties, not methods
    x = embedding.x
    y = embedding.y
    s = embedding.s
    t = embedding.t
    
    # Determine molecule groups (trees from the same original molecule)
    molecule_ids = []
    for idx, tree in enumerate(trees):
        # Extract molecule ID from the root label (expected format: "tree{mol_idx}_{tree_idx}")
        label = tree.root.label
        if "_" in label:
            mol_id = label.split("_")[0].replace("tree", "")
        else:
            mol_id = "0"  # Default if label doesn't follow expected format
        molecule_ids.append(mol_id)
    
    unique_mol_ids = list(set(molecule_ids))
    mol_id_to_idx = {mol_id: idx for idx, mol_id in enumerate(unique_mol_ids)}
    mol_indices = [mol_id_to_idx[mol_id] for mol_id in molecule_ids]
    
    # Plot the layout
    print("Plotting layout...")
    plt.figure(figsize=(12, 10))
    
    # Create a custom colormap with good distinguishable colors
    cmap = plt.get_cmap('tab10', max(10, len(unique_mol_ids)))
    
    # Plot edges
    for i in range(len(s)):
        plt.plot(
            [x[s[i]], x[t[i]]],
            [y[s[i]], y[t[i]]],
            "k-",
            linewidth=0.5,
            alpha=0.3,
            zorder=1,
        )
    
    # Plot vertices with size and color based on properties
    # Use molecule ID for color
    plt.scatter(
        x, y, 
        s=80,  # Fixed size for all nodes
        c=mol_indices, 
        cmap=cmap, 
        alpha=0.8,
        edgecolors='black',
        linewidths=0.5,
        zorder=2
    )
    
    # Add labels for each point (using tree index within molecule)
    for i, (mol_id, tree) in enumerate(zip(molecule_ids, trees)):
        # Extract tree index from the root label
        label = tree.root.label
        if "_" in label:
            tree_idx = label.split("_")[1]
        else:
            tree_idx = "0"  # Default if label doesn't follow expected format
            
        plt.annotate(
            f"{tree_idx}",
            (x[i], y[i]),
            textcoords="offset points",
            xytext=(0, 0),
            ha="center",
            va="center",
            fontsize=8,
            color="white",
            fontweight="bold",
        )
    
    # Add a legend mapping colors to molecule IDs
    handles = []
    labels = []
    for mol_id in unique_mol_ids:
        idx = mol_id_to_idx[mol_id]
        patch = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=cmap(idx), 
                          markersize=10, label=f"Mol {mol_id}")
        handles.append(patch)
        labels.append(f"Mol {mol_id}")
    
    plt.legend(handles, labels, loc='upper right', title="Molecules")
    
    plt.title("Synthesis Route Trees Grouped by Tree Edit Distance")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(args.output, dpi=300)
    print(f"Figure saved as '{args.output}'")

if __name__ == "__main__":
    main() 