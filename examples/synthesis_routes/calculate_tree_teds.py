#!/usr/bin/env python
"""
Script to calculate TED (Tree Edit Distance) between specific trees from an AiZynthFinder output file.
This uses the exact same calculation methods as the tmap visualization code.
"""

import os
import sys
import argparse
import gzip
import json
import pandas as pd
import numpy as np
from typing import List, Dict, Any, Tuple, Optional, Union

import use_local_tmap
from tmap.core.synthesis_tree import SynthesisTree, SynthesisTreeNode
from tmap.helpers.ted_utils import compute_ted, convert_to_apted_tree, compute_ted_with_apted
from tmap.helpers.ted_utils import get_reaction_classifications, CustomConfig

# Check if apted is available
try:
    from apted import APTED, Config
    APTED_AVAILABLE = True
except ImportError:
    APTED_AVAILABLE = False
    print("Warning: APTED package not found. Install with 'pip install apted'")


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

def load_aizynthfinder_data(filepath):
    """
    Load trees from an AiZynthFinder output file.
    
    Args:
        filepath: Path to the AiZynthFinder output.json.gz file
        
    Returns:
        pd.DataFrame: DataFrame with AiZynthFinder data
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
            return None
    
    return data

def print_tree_structure(tree_data, level=0, max_level=3):
    """Print a simplified view of the tree structure for debugging"""
    if level > max_level:
        print("  " * level + "... (truncated)")
        return
    
    # Print current node
    node_type = tree_data.get("type", "unknown")
    if node_type == "mol":
        smiles = tree_data.get("smiles", "")
        label = f"Molecule: {smiles[:20]}..." if len(smiles) > 20 else f"Molecule: {smiles}"
    elif node_type == "reaction":
        classification = tree_data.get("metadata", {}).get("classification", "")
        label = f"Reaction: {classification}"
    else:
        label = f"Node: {node_type}"
    
    print("  " * level + label)
    
    # Print children
    children = tree_data.get("children", [])
    for child in children:
        print_tree_structure(child, level + 1, max_level)

def main():
    parser = argparse.ArgumentParser(description="Calculate TED between specific trees in AiZynthFinder output")
    parser.add_argument("input_file", help="Path to the AiZynthFinder output.json.gz file")
    parser.add_argument("--mol1", type=int, default=0, help="Index of first molecule (default: 0)")
    parser.add_argument("--tree1", type=int, default=0, help="Index of first tree within molecule (default: 0)")
    parser.add_argument("--mol2", type=int, default=0, help="Index of second molecule (default: 0)")
    parser.add_argument("--tree2", type=int, default=1, help="Index of second tree within molecule (default: 1)")
    parser.add_argument("--print-trees", action="store_true", help="Print simplified tree structures")
    parser.add_argument("--verbose", action="store_true", help="Print detailed information during TED calculation")
    parser.add_argument("--matrix", action="store_true", help="Calculate TED matrix for all trees of the specified molecules")
    
    args = parser.parse_args()
    
    # Check if the input file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        return
    
    # Load the trees from the AiZynthFinder output file
    print(f"Loading trees from {args.input_file}...")
    data = load_aizynthfinder_data(args.input_file)
    
    if data is None or not hasattr(data, "trees"):
        print("No trees found in the input file.")
        return
    
    # Calculate TED matrix for all trees of specified molecules
    if args.matrix:
        if len(data.trees) <= args.mol1:
            print(f"Error: Molecule index {args.mol1} is out of range. Only {len(data.trees)} molecules available.")
            return
            
        trees1 = data.trees[args.mol1]
        print(f"Found {len(trees1)} trees for molecule {args.mol1}")
        
        # Convert trees to SynthesisTree objects
        synthesis_trees = []
        for i, tree_data in enumerate(trees1):
            synthesis_tree = convert_aizynthfinder_tree_to_synthesis_tree(tree_data, tree_idx=f"{args.mol1}_{i}")
            synthesis_trees.append(synthesis_tree)
        
        # Calculate TED matrix
        n = len(synthesis_trees)
        ted_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                ted = compute_ted(synthesis_trees[i], synthesis_trees[j])
                ted_matrix[i, j] = ted
                ted_matrix[j, i] = ted
        
        # Print TED matrix
        print("\nTED Matrix:")
        print("    " + " ".join(f"{i:5d}" for i in range(n)))
        for i in range(n):
            print(f"{i:3d} " + " ".join(f"{ted_matrix[i, j]:5.2f}" for j in range(n)))
        
        # Print statistics
        print("\nTED Statistics:")
        triu_indices = np.triu_indices(n, k=1)  # Upper triangular indices, excluding diagonal
        ted_values = ted_matrix[triu_indices]
        print(f"Min TED: {np.min(ted_values):.2f}")
        print(f"Max TED: {np.max(ted_values):.2f}")
        print(f"Mean TED: {np.mean(ted_values):.2f}")
        print(f"Median TED: {np.median(ted_values):.2f}")
        
        return
    
    # Check if the molecule indices are valid
    if len(data.trees) <= args.mol1:
        print(f"Error: First molecule index {args.mol1} is out of range. Only {len(data.trees)} molecules available.")
        return
    
    if len(data.trees) <= args.mol2:
        print(f"Error: Second molecule index {args.mol2} is out of range. Only {len(data.trees)} molecules available.")
        return
    
    # Get the trees for the specified molecules
    trees1 = data.trees[args.mol1]
    trees2 = data.trees[args.mol2]
    
    # Check if the tree indices are valid
    if len(trees1) <= args.tree1:
        print(f"Error: First tree index {args.tree1} is out of range. Only {len(trees1)} trees available for molecule {args.mol1}.")
        return
    
    if len(trees2) <= args.tree2:
        print(f"Error: Second tree index {args.tree2} is out of range. Only {len(trees2)} trees available for molecule {args.mol2}.")
        return
    
    # Get the specific trees
    tree1_data = trees1[args.tree1]
    tree2_data = trees2[args.tree2]
    
    # Print tree structures if requested
    if args.print_trees:
        print(f"\nTree structure for molecule {args.mol1}, tree {args.tree1}:")
        print_tree_structure(tree1_data)
        
        print(f"\nTree structure for molecule {args.mol2}, tree {args.tree2}:")
        print_tree_structure(tree2_data)
    
    # Convert to SynthesisTree objects
    print(f"\nConverting trees to SynthesisTree objects...")
    synthesis_tree1 = convert_aizynthfinder_tree_to_synthesis_tree(tree1_data, tree_idx=f"{args.mol1}_{args.tree1}")
    synthesis_tree2 = convert_aizynthfinder_tree_to_synthesis_tree(tree2_data, tree_idx=f"{args.mol2}_{args.tree2}")
    
    # Enable verbose debugging for APTED
    if args.verbose:
        # Convert to APTED format and print
        apted_tree1 = convert_to_apted_tree(synthesis_tree1)
        apted_tree2 = convert_to_apted_tree(synthesis_tree2)
        
        print("\nAPTED Tree 1:")
        print(json.dumps(apted_tree1, indent=2))
        
        print("\nAPTED Tree 2:")
        print(json.dumps(apted_tree2, indent=2))
        
        # Create a debugging version of CustomConfig
        class DebugCustomConfig(CustomConfig):
            def rename(self, node1, node2):
                cost = super().rename(node1, node2)
                print(f"Rename: {node1.get('type')} -> {node2.get('type')}, Cost: {cost}")
                if node1.get('type') == 'mol':
                    smiles1 = node1.get('metadata', {}).get('smiles', node1.get('label', ''))
                    smiles2 = node2.get('metadata', {}).get('smiles', node2.get('label', ''))
                    print(f"  Mol1: {smiles1}")
                    print(f"  Mol2: {smiles2}")
                if node1.get('type') == 'reaction':
                    class1 = node1.get('metadata', {}).get('classification', '')
                    class2 = node2.get('metadata', {}).get('classification', '')
                    print(f"  Reaction1: {class1}")
                    print(f"  Reaction2: {class2}")
                return cost
                
            def delete(self, node):
                cost = super().delete(node)
                print(f"Delete: {node.get('type')}, Cost: {cost}")
                return cost
                
            def insert(self, node):
                cost = super().insert(node)
                print(f"Insert: {node.get('type')}, Cost: {cost}")
                return cost
        
        # Calculate TED with debugging
        if APTED_AVAILABLE:
            print("\nCalculating TED with verbose debugging...")
            apted = APTED(apted_tree1, apted_tree2, DebugCustomConfig())
            ted = apted.compute_edit_distance()
            print(f"\nTED: {ted}")
        else:
            print("\nAPTED not available, cannot calculate TED with verbose debugging")
    
    # Calculate TED
    print("\nCalculating TED...")
    ted = compute_ted(synthesis_tree1, synthesis_tree2)
    print(f"\nTree Edit Distance (TED): {ted}")
    
    return ted

if __name__ == "__main__":
    main() 