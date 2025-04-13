"""
Simple example using tmap to visualize synthesis trees
Following the pattern of the MNIST example from tmap documentation
"""
import os
import sys
import numpy as np
import tmap as tm
from faerun import Faerun

# Add module paths (ensure tmap modules can be found)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT_DIR = os.path.abspath(os.path.join(_SCRIPT_DIR, '../..'))
sys.path.insert(0, os.path.join(_ROOT_DIR, 'src'))

# Create some simple data - 10 groups of points in 2D space
np.random.seed(42)
n_groups = 10
n_points_per_group = 20
n_points = n_groups * n_points_per_group

# Create the data points
data = []
labels = []
group_ids = []

for i in range(n_groups):
    # Create a cluster of points
    center_x = np.random.uniform(-10, 10)
    center_y = np.random.uniform(-10, 10)
    
    for j in range(n_points_per_group):
        # Add some noise to create a cluster
        x = center_x + np.random.normal(0, 1)
        y = center_y + np.random.normal(0, 1)
        
        # Store as a feature vector
        data.append([x, y])
        labels.append(f"Group {i}, Point {j}")
        group_ids.append(i)

# Convert to numpy arrays
data = np.array(data)
group_ids = np.array(group_ids)

# Create a configuration for the layout
cfg = tm.LayoutConfiguration()
cfg.k = 10                      # Number of nearest neighbors
cfg.kc = 10                     # Number to multiply k by for nearest neighbor search
cfg.fme_iterations = 500        # Number of iterations of the force-directed layout
cfg.sl_repeats = 2              # Number of repeats of the scaling layout algorithm

# Create a weight matrix (distances between points)
n = len(data)
weights = np.zeros((n, n), dtype=np.float32)

# Compute the Euclidean distances between points
for i in range(n):
    for j in range(i+1, n):
        dist = np.linalg.norm(data[i] - data[j])
        weights[i, j] = dist
        weights[j, i] = dist

# Create an edge list from the weight matrix
edges = []
for i in range(n):
    # Find the k nearest neighbors for each point
    nn_indices = np.argsort(weights[i])[:cfg.k+1]  # +1 because first one is self
    for j in nn_indices:
        if i != j:  # Skip self
            edges.append((i, j, weights[i, j]))

# Generate the layout using tmap
x, y, s, t, _ = tm.layout_from_edge_list(n, edges, cfg)

# Convert the lists to NumPy arrays for Faerun
x = np.array(x)
y = np.array(y)

# Create the Faerun visualization
f = Faerun(clear_color="#222222", view="front", coords=False)

# Add a scatter plot
f.add_scatter(
    "simple_scatter",
    {
        "x": x,
        "y": y,
        "c": group_ids,         # Color by group ID (integer)
        "labels": labels,       # Add labels for tooltips
    },
    shader="smoothCircle",      # Use a smooth circle shader
    point_scale=5.0,            # Size of the points
    categorical=True,           # The colors are categorical
    colormap="tab10",           # Use tab10 colormap
    has_legend=True,            # Add a legend
    legend_labels=[(i, f"Group {i}") for i in range(n_groups)],  # Legend entries
)

# Add edges
edges_data = {
    "from": s,
    "to": t
}

# Add the edges as a tree
f.add_tree(
    "simple_tree",
    edges_data,
    point_helper="simple_scatter",  # Use the points from the scatter plot
    color="#666666",               # Gray edges
)

# Save the visualization
output_path = "simple_tmap_viz"
f.plot(output_path)
print(f"Visualization saved to {output_path}.html") 