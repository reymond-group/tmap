"""
Example of using the TEDLayoutGenerator to visualize synthesis route trees from AiZynthFinder output.
This script loads trees from an AiZynthFinder output.json.gz file and visualizes them
using Tree Edit Distance (TED) to group similar synthesis routes together.
Uses Faerun for interactive web-based visualization.
"""
import use_local_tmap
import numpy as np
import pandas as pd
from faerun import Faerun
import argparse
import os
import json
import gzip
from matplotlib.colors import ListedColormap
from tmap.core.synthesis_tree import SynthesisTree, SynthesisTreeNode
from tmap.layout_generators.ted_layout_generator import TEDLayoutGenerator
from tmap.helpers.ted_utils import compute_ted
import colorcet as cc

colors = cc.glasbey


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
        smiles = node_data.get("smiles", "")
        # Store molecule ID if available in metadata
        mol_id = node_data.get("metadata", {}).get("molecule_id", "")
        label = f"mol-{mol_id}-{parent_label}" if mol_id else smiles if smiles else f"Molecule-{parent_label}"
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
            "smiles": node_data.get("smiles", ""),
            "mol_id": node_data.get("metadata", {}).get("molecule_id", ""),
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

def extract_smiles_from_tree(tree):
    """
    Extract SMILES from a tree if available.
    
    Args:
        tree: A SynthesisTree object
        
    Returns:
        str: SMILES string if available, otherwise empty string
    """
    label = tree.root.label
    if label.startswith("Reaction"):
        return ""
    
    # For molecule nodes, the label might be a SMILES string
    if "." not in label and "-" not in label:
        return label
    return ""

def main():
    """
    Main function to demonstrate the TEDLayoutGenerator with AiZynthFinder trees using Faerun.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Visualize AiZynthFinder trees using Tree Edit Distance and Faerun")
    parser.add_argument("input_file", help="Path to the AiZynthFinder output.json.gz file")
    parser.add_argument("--k", type=int, default=5, help="Number of nearest neighbors (default: 5)")
    parser.add_argument("--ted-threshold", type=float, default=5.0, 
                      help="Maximum Tree Edit Distance (TED) allowed between connected nodes. Lower values = more strict similarity (default: 5.0)")
    parser.add_argument("--force-connect", action="store_true",
                      help="Force connection between all disconnected components, may break TED threshold")
    parser.add_argument("--title", default="Synthesis Route Trees", help="Title for the visualization (default: Synthesis Route Trees)")
    parser.add_argument("--output", default="aizynthfinder_trees_faerun", help="Output HTML file name without extension (default: aizynthfinder_trees_faerun)")
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
    print(f"TED threshold: {args.ted_threshold}")
    print(f"Force connect components: {args.force_connect}")
    
    layout_generator = TEDLayoutGenerator(
        k=args.k,                      # Number of nearest neighbors
        create_mst=True,               # Create minimum spanning tree
        keep_knn=True,                 # Keep k-nearest neighbors graph
        ted_threshold=args.ted_threshold,     # TED threshold
        force_connect=args.force_connect      # Force connection between components
    )
    
    # Generate the layout
    print("Generating layout...")
    embedding = layout_generator.layout(trees)
    
    # Extract coordinates and edges from embedding properties
    x = embedding.x
    y = embedding.y
    s = embedding.s
    t = embedding.t
    
    # Convert the edges to the format expected by Faerun
    edges = []
    for i in range(len(s)):
        edges.append([int(s[i]), int(t[i])])
    
    # Determine molecule groups (trees from the same original molecule)
    molecule_ids = []
    labels = []
    smiles = []
    node_types = []
    
    # Map to track which index corresponds to which molecule
    mol_id_mapping = {}
    next_mol_id = 0
    
    print(f"Analyzing {len(trees)} trees to extract molecule data...")
    
    for idx, tree in enumerate(trees):
        # Try to extract molecule data from the tree
        root = tree.root
        root_data = root.data
        root_type = root_data.get("type", "unknown")
        root_smiles = root_data.get("smiles", "")
        
        # Check if we have molecule ID in the data
        mol_id_from_data = root_data.get("mol_id", "")
        
        # Use SMILES as a unique identifier if available, otherwise use the label
        if root_smiles and len(root_smiles) > 0:
            mol_key = root_smiles
        else:
            # Try to extract from label
            label = root.label
            if label.startswith("mol-") and "-" in label:
                parts = label.split("-")
                if len(parts) > 1 and parts[1].isdigit():
                    mol_id_from_data = parts[1]
                mol_key = f"molecule-{mol_id_from_data}" if mol_id_from_data else label
            else:
                mol_key = label
        
        # Assign a unique ID to each molecule
        if mol_key not in mol_id_mapping:
            mol_id_mapping[mol_key] = next_mol_id
            next_mol_id += 1
        
        # Get the numeric ID for this molecule
        mol_numeric_id = mol_id_mapping[mol_key]
        molecule_ids.append(mol_numeric_id)
        
        # Create a descriptive label for the node
        if mol_id_from_data:
            tree_label = f"Mol {mol_id_from_data}, Tree {idx}"
        else:
            tree_label = f"Mol {mol_numeric_id}, Tree {idx}"
        labels.append(tree_label)
        
        # Store other data
        smiles.append(root_smiles)
        node_types.append(root_type)
    
    print(f"Found {len(mol_id_mapping)} unique molecules.")
    
    # Create color data based on molecule ID
    # Convert to numpy array to ensure all values are numerical
    color_data = np.array(molecule_ids, dtype=np.int32)
    
    # Create category data based on node type
    category_data = ["Molecule" if t == "mol" else "Reaction" for t in node_types]
    
    # Create Faerun visualization
    print("Creating Faerun visualization...")
    f = Faerun(clear_color="#222222", view="front", coords=False, title=args.title)
    
    # Get unique molecule IDs for legend labels
    unique_mol_ids = sorted(set(molecule_ids))
    legend_labels = [(i, f"Molecule {i}") for i in unique_mol_ids]
    
    # Create a new plot
    f.add_scatter(
        "aizynthfinder_trees",
        {
            "x": x,
            "y": y,
            "c": color_data,
            "labels": labels,
        },
        shader="smoothCircle",
        colormap=ListedColormap(colors),
        point_scale=5.0,
        has_legend=True,
        legend_labels=legend_labels,
        categorical=True,
        interactive=True
    )
    
    # Add edges
    f.add_tree(
        "tree",
        {"from": [e[0] for e in edges], "to": [e[1] for e in edges]},
        point_helper="aizynthfinder_trees",
        color="#666666",
    )
    
    # Save the visualization to HTML
    print(f"Saving visualization to {args.output}.html...")
    f.plot(args.output)
    print(f"Visualization saved as {args.output}.html")

if __name__ == "__main__":
    main() 