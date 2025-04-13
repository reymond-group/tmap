# Synthesis Route Visualization with TMAP

This folder contains tools for visualizing chemical synthesis routes using TMAP (Tree-based Map) and Tree Edit Distance (TED) to group similar synthesis trees together.

## Overview

The tools in this folder enable you to:

1. Load synthesis route trees from AiZynthFinder output files
2. Group similar trees based on Tree Edit Distance (TED)
3. Visualize the trees in a 2D layout where similar trees are positioned close to each other

## Setup

Install the required packages:

```bash
pip install -r requirements.txt
```

## Usage

### Visualizing AiZynthFinder Output

To visualize synthesis routes from an AiZynthFinder output file:

```bash
python aizynthfinder_visualization.py path/to/output.json.gz --k 5 --output visualization.png
```

Parameters:
- `path/to/output.json.gz`: Path to the AiZynthFinder output file
- `--k`: Number of nearest neighbors to use (default: 5)
- `--output`: Output image file path (default: aizynthfinder_trees_ted.png)

### Example with Sample Trees

To run the simple example with generated sample trees:

```bash
python synthesis_tree_example.py
```

## Understanding the Visualization

In the visualization:
- Each point represents a synthesis route tree
- Points with the same color belong to the same molecule
- Numbers inside points indicate the tree index within that molecule
- Lines connect similar trees (based on Tree Edit Distance)
- Trees that are structurally similar are positioned closer together

## Tree Edit Distance

Tree Edit Distance (TED) measures the similarity between two trees as the minimum cost of transforming one tree into another using operations like node deletion, insertion, and renaming.

The TED calculation uses custom costs for synthetic routes:
- Different node types (mol vs. reaction) are considered highly different
- For reaction nodes, differences in reaction classification are weighted:
  - First-level classification difference: 0.8
  - Second-level classification difference: 0.5
  - Third-level classification difference: 0
- Molecule nodes are considered similar (renaming cost = 0)

## Technical Details

The visualization uses:
- `TEDLayoutGenerator`: A custom TMAP layout generator that uses TED for measuring tree similarity
- `CustomConfig`: A custom APTED configuration for computing Tree Edit Distance with reaction-specific metrics
- OGDF (Open Graph Drawing Framework) for the actual graph layout

The TMAP algorithm creates a k-nearest neighbor graph based on TED distances, computes a minimum spanning tree, and applies a force-directed layout to position the nodes in 2D space. 