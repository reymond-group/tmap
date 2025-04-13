from typing import List, Dict, Union
import json
import subprocess
import os

def extract_pathway_graph(node: Dict, parent_smiles=None):
    """Extract a graph representation of the synthesis pathway.
    
    Returns a tuple of (nodes, edges) where:
    - nodes is a dict mapping node_id to node properties
    - edges is a list of (source_id, target_id) tuples
    """
    nodes = {}
    edges = []
    
    # Process current node
    if node.get("type") == "mol":
        # This is a molecule node
        mol_id = f"mol_{id(node)}"
        mol_smiles = node.get("smiles", "")
        in_stock = node.get("in_stock", False)
        
        nodes[mol_id] = {
            "type": "molecule",
            "smiles": mol_smiles,
            "in_stock": in_stock
        }
        
        # Connect to parent if it exists
        if parent_smiles:
            edges.append((parent_smiles, mol_id))
    
    elif node.get("type") == "reaction":
        # This is a reaction node
        reaction_id = f"rxn_{id(node)}"
        reaction_smiles = node.get("smiles", "")
        
        # Get reaction name/classification from metadata if available
        reaction_name = "Unknown Reaction"
        if "metadata" in node and "classification" in node["metadata"]:
            reaction_name = node["metadata"]["classification"]
            
        nodes[reaction_id] = {
            "type": "reaction",
            "name": reaction_name,
            "smiles": reaction_smiles
        }
        
        # Connect to parent if it exists
        if parent_smiles:
            edges.append((parent_smiles, reaction_id))
        
        # Set as parent for children
        parent_smiles = reaction_id
    
    # Process children recursively
    for child in node.get("children", []):
        child_nodes, child_edges = extract_pathway_graph(child, parent_smiles)
        nodes.update(child_nodes)
        edges.extend(child_edges)
    
    return nodes, edges

def generate_graphviz_dot(nodes, edges):
    """Generate a Graphviz DOT representation of the synthesis pathway."""
    dot = ["digraph SynthesisPathway {"]
    
    # Set global graph attributes
    dot.append("  rankdir=LR;")  # Left-to-right layout
    dot.append("  node [shape=box, style=filled, fontname=\"Arial\"];")
    
    # Add nodes
    for node_id, props in nodes.items():
        if props["type"] == "molecule":
            label = props["smiles"]
            if len(label) > 30:
                label = label[:27] + "..."
            color = "#C6EFCE" if props["in_stock"] else "#FFEB9C"  # Green if in stock, yellow if not
            dot.append(f"  \"{node_id}\" [label=\"{label}\", fillcolor=\"{color}\"];")
        else:  # reaction
            label = props["name"]
            dot.append(f"  \"{node_id}\" [label=\"{label}\", shape=ellipse, fillcolor=\"#B4C6E7\"];")  # Blue for reactions
    
    # Add edges
    for source, target in edges:
        dot.append(f"  \"{source}\" -> \"{target}\";")
    
    dot.append("}")
    return "\n".join(dot)

def save_pathway_representation(root: Dict, output_dot="synthesis_pathway.dot", output_json="synthesis_pathway.json", output_png="synthesis_pathway.png"):
    """Extract the pathway and save it in DOT, JSON, and PNG formats."""
    
    # Extract the pathway
    nodes, edges = extract_pathway_graph(root)
    
    # Save as DOT format for visualization
    dot_content = generate_graphviz_dot(nodes, edges)
    with open(output_dot, "w") as f:
        f.write(dot_content)
    
    # Save as JSON for further processing
    pathway_data = {
        "nodes": nodes,
        "edges": edges
    }
    with open(output_json, "w") as f:
        json.dump(pathway_data, f, indent=2)
    
    # Generate PNG using Graphviz's dot command
    generate_png_from_dot(output_dot, output_png)
    
    return dot_content, pathway_data

def generate_png_from_dot(dot_file, png_file):
    """Generate a PNG image from a DOT file using Graphviz."""
    try:
        # Run the dot command to generate the PNG
        subprocess.run(["dot", "-Tpng", dot_file, "-o", png_file], check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error generating PNG: {e}")
        return False
    except FileNotFoundError:
        print("Graphviz's 'dot' command not found. Please install Graphviz to generate PNG images.")
        print("Installation: https://graphviz.org/download/")
        return False

# Example usage
# Load the reaction tree from JSON file
with open("test.json") as f:
    reaction_tree = json.load(f)

dot_content, pathway_data = save_pathway_representation(reaction_tree)
print(f"Synthesis pathway extracted:")
print(f"- {len(pathway_data['nodes'])} nodes")
print(f"- {len(pathway_data['edges'])} edges")
print(f"DOT file saved to 'synthesis_pathway.dot'")
print(f"JSON representation saved to 'synthesis_pathway.json'")

# Check if PNG was successfully generated
if os.path.exists("synthesis_pathway.png"):
    print(f"PNG image saved to 'synthesis_pathway.png'")
else:
    print("\nTo manually visualize the pathway using Graphviz, run:")
    print("$ dot -Tpng synthesis_pathway.dot -o synthesis_pathway.png")
