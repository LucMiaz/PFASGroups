# generate graph of HalogenGroups based on file test/specificity_test_groups.json
import json
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np
import os
current_dir = os.path.dirname(os.path.abspath(__file__))

def load_graph_from_json(filename):
    """Load the list-based JSON structure and convert to NetworkX directed graph with edge types."""
    with open(filename, 'r') as f:
        data = json.load(f)
    
    G = nx.DiGraph()
    
    # Process each edge in the list
    for edge_data in data:
        source = edge_data['source']
        target = edge_data['target']
        edge_type = edge_data['edge_type']
        
        # Add nodes if they don't exist
        G.add_node(source)
        G.add_node(target)
        
        # Add edge with type information
        G.add_edge(source, target, edge_type="oneway")
        if edge_type == 'poly':
            G.add_edge(target, source, edge_type='poly')
        elif edge_type == 'per':
            G.add_edge(target, source, edge_type='per')
    return G

def visualize_graph(G, title="PFAS Groups Network"):
    """Create a visualization of the graph with hierarchical layout and edge types."""
    plt.figure(figsize=(20, 16))
    
    # Use hierarchical layout (graphviz_layout requires pygraphviz)
    try:
        pos = nx.nx_agraph.graphviz_layout(G, prog='neato')
    except:
        # Fallback to spring layout if graphviz is not available
        print("Graphviz not available, using spring layout...")
        pos = nx.spring_layout(G, k=4, iterations=50, seed=40)
    
    # Separate edges by type for different styling
    per_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'per']
    poly_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'poly']
    oneway_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'oneway']
    
    # Draw edges with different colors and styles based on type
    if per_edges:
        nx.draw_networkx_edges(G, pos, edgelist=per_edges, edge_color='#01665e', 
                      arrows=True, arrowsize=15, arrowstyle='->', 
                      alpha=0.7, width=2, style='dashed', connectionstyle='arc3,rad=0.1',
                      label='Co-occurs if perfluorinated')
    
    if poly_edges:
        nx.draw_networkx_edges(G, pos, edgelist=poly_edges, edge_color='#8c510a', 
                              arrows=True, arrowsize=15, arrowstyle='->', connectionstyle='arc3,rad=0.1',
                              alpha=0.7, width=2, style='dashed', label='Co-occurs if polyfluorinated')
    if oneway_edges:
        nx.draw_networkx_edges(G, pos, edgelist=poly_edges, edge_color='#8c510a', 
                              arrows=True, arrowsize=15, arrowstyle='->', connectionstyle='arc3,rad=0.1',
                              alpha=0.7, width=2, style='solid', label='Co-occurs with')
    
    # Categorize nodes by their characteristics for coloring
    node_colors = []
    node_sizes = []
    for node in G.nodes():
        # Color nodes based on their names
        if 'Perfluoro' in node:
            node_colors.append('#b3e2cd')  # Red for perfluoro compounds
            node_sizes.append(800)
        elif 'Polyfluoro' in node:
            node_colors.append('#fdcdac')  # Teal for polyfluoro compounds
            node_sizes.append(700)
        elif 'acid' in node.lower():
            node_colors.append('#cbd5e8')  # Blue for acids
            node_sizes.append(600)
        elif 'alcohol' in node.lower():
            node_colors.append('#f4cae4')  # Green for alcohols
            node_sizes.append(600)
        elif 'ether' in node.lower():
            node_colors.append('#e6f5c9')  # Yellow for ethers
            node_sizes.append(600)
        else:
            node_colors.append('#fff2ae')  # Purple for others
            node_sizes.append(500)
    
    # Draw nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                          alpha=0.8, edgecolors='black', linewidths=1)
    
    # Draw labels with better formatting
    labels = {}
    for node in G.nodes():
        # Wrap long labels
        if len(node) > 20:
            words = node.split(' ')
            if len(words) > 1:
                mid = len(words) // 2
                labels[node] = ' '.join(words[:mid]) + '\n' + ' '.join(words[mid:])
            else:
                labels[node] = node[:20] + '\n' + node[20:]
        else:
            labels[node] = node
    
    nx.draw_networkx_labels(G, pos, labels, font_size=12)
    
    plt.title(title, fontsize=16, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    
    # Add comprehensive legend
    legend_elements = [
        # Node types
        plt.scatter([], [], c='#b3e2cd', s=100, label='Perfluoro compounds', marker='o'),
        plt.scatter([], [], c='#fdcdac', s=100, label='Polyfluoro compounds', marker='o'),
        plt.scatter([], [], c='#cbd5e8', s=100, label='Acids', marker='o'),
        plt.scatter([], [], c='#f4cae4', s=100, label='Alcohols', marker='o'),
        plt.scatter([], [], c='#e6f5c9', s=100, label='Ethers', marker='o'),
        plt.scatter([], [], c='#fff2ae', s=100, label='Other compounds', marker='o'),
        # Edge types
        plt.Line2D([0], [0], color='#8c510a', linewidth=2, linestyle='-', label='Co-occurs'),
        plt.Line2D([0], [0], color='#01665e', linewidth=2, linestyle='--', label='Co-occurs if perfluorinated'),
        plt.Line2D([0], [0], color='#8c510a', linewidth=2, linestyle='--', label='Co-occurs if polyfluorinated'),
    ]
    plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1))
    
    return plt

def analyze_graph(G):
    """Print basic statistics about the graph including edge types."""
    print(f"Graph Statistics:")
    print(f"- Number of nodes: {G.number_of_nodes()}")
    print(f"- Number of edges: {G.number_of_edges()}")
    print(f"- Is connected: {nx.is_weakly_connected(G)}")
    print(f"- Number of weakly connected components: {nx.number_weakly_connected_components(G)}")
    
    # Analyze edge types
    edge_types = {}
    for u, v, d in G.edges(data=True):
        edge_type = d.get('edge_type', 'unknown')
        edge_types[edge_type] = edge_types.get(edge_type, 0) + 1
    
    print(f"- Edge type distribution:")
    for edge_type, count in edge_types.items():
        print(f"  - {edge_type}: {count} edges")
    
    # Find root nodes (nodes with no predecessors)
    root_nodes = [node for node in G.nodes() if G.in_degree(node) == 0]
    print(f"- Root nodes ({len(root_nodes)}): {root_nodes}")
    
    # Find leaf nodes (nodes with no successors)
    leaf_nodes = [node for node in G.nodes() if G.out_degree(node) == 0]
    print(f"- Leaf nodes ({len(leaf_nodes)}): {leaf_nodes}")
    
    # Find most connected nodes
    if G.nodes():
        max_in_degree = max(G.in_degree(node) for node in G.nodes())
        max_out_degree = max(G.out_degree(node) for node in G.nodes())
        
        most_referenced = [node for node in G.nodes() if G.in_degree(node) == max_in_degree]
        most_connecting = [node for node in G.nodes() if G.out_degree(node) == max_out_degree]
        
        print(f"- Most referenced nodes (in-degree {max_in_degree}): {most_referenced}")
        print(f"- Most connecting nodes (out-degree {max_out_degree}): {most_connecting}")
    
    # Show some example paths with edge types
    print(f"\n- Sample paths with edge types:")
    
    # Find some interesting paths
    paths_shown = 0
    max_paths = 10
    
    for source in list(G.nodes())[:5]:  # Try first 5 nodes as sources
        for target in list(G.nodes()):
            if source != target and nx.has_path(G, source, target):
                try:
                    path = nx.shortest_path(G, source, target)
                    if len(path) >= 3:  # Only show paths with at least 2 edges
                        path_description = []
                        for i in range(len(path) - 1):
                            edge_type = G[path[i]][path[i+1]].get('edge_type', 'unknown')
                            path_description.append(f"{path[i]} --({edge_type})--> {path[i+1]}")
                        
                        print(f"  {' → '.join(path_description)}")
                        paths_shown += 1
                        
                        if paths_shown >= max_paths:
                            break
                except nx.NetworkXNoPath:
                    continue
        if paths_shown >= max_paths:
            break
    
    # Show nodes that appear in both per and poly relationships
    print(f"\n- Nodes with both perfluorinated and polyfluorinated relationships:")
    nodes_with_both = {}
    for u, v, d in G.edges(data=True):
        edge_type = d.get('edge_type', 'unknown')
        if v not in nodes_with_both:
            nodes_with_both[v] = set()
        nodes_with_both[v].add(edge_type)
    
    dual_nodes = [(node, types) for node, types in nodes_with_both.items() if 'per' in types and 'poly' in types]
    for node, types in dual_nodes[:10]:  # Show first 10
        print(f"  - {node}: {sorted(types)}")
    
    if len(dual_nodes) > 10:
        print(f"  ... and {len(dual_nodes) - 10} more")
    elif len(dual_nodes) == 0:
        print(f"  No nodes found with both edge types")

def analyze_edge_patterns(G):
    """Analyze patterns in the edge relationships."""
    print(f"\n=== Edge Pattern Analysis ===")
    
    # Count parallel edges (same source-target pair with different edge types)
    edge_pairs = {}
    for u, v, d in G.edges(data=True):
        pair = (u, v)
        if pair not in edge_pairs:
            edge_pairs[pair] = []
        edge_pairs[pair].append(d.get('edge_type', 'unknown'))
    
    parallel_edges = [(pair, types) for pair, types in edge_pairs.items() if len(types) > 1]
    print(f"- Parallel edges (same nodes, different types): {len(parallel_edges)}")
    
    if parallel_edges:
        print("  Examples:")
        for (source, target), types in parallel_edges[:5]:
            print(f"    {source} → {target}: {sorted(set(types))}")
    
    # Analyze perfluorinated vs polyfluorinated preference by node
    per_nodes = set()
    poly_nodes = set()
    
    for u, v, d in G.edges(data=True):
        edge_type = d.get('edge_type', 'unknown')
        if edge_type == 'per':
            per_nodes.add(u)
            per_nodes.add(v)
        elif edge_type == 'poly':
            poly_nodes.add(u)
            poly_nodes.add(v)
    
    per_only = per_nodes - poly_nodes
    poly_only = poly_nodes - per_nodes
    both_types = per_nodes & poly_nodes
    
    print(f"- Nodes in perfluorinated pathways only: {len(per_only)}")
    print(f"- Nodes in polyfluorinated pathways only: {len(poly_only)}")
    print(f"- Nodes in both pathway types: {len(both_types)}")

# Main execution
filename = f'{current_dir}/specificity_test_groups.json'

if __name__ == "__main__":
    try:
        # Load and create the graph
        print(f"Loading graph from {filename}...")
        G = load_graph_from_json(filename)
        
        # Analyze the graph
        analyze_graph(G)
        analyze_edge_patterns(G)
        print()
        
        # Visualize the graph
        print("Creating visualization...")
        plt = visualize_graph(G, "PFAS Chemical Groups Hierarchy")
        
        # Save the plot
        output_filename = filename.replace('.json', '_network.pdf')
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Graph saved as {output_filename}")
        
        # Show the plot
        plt.show()
        
        # Optional: Save the graph in other formats
        print("\nSaving graph in additional formats...")
        
        # Save as GraphML (can be opened in Gephi, Cytoscape, etc.)
        graphml_filename = filename.replace('.json', '_network.graphml')
        nx.write_graphml(G, graphml_filename)
        print(f"Graph saved as GraphML: {graphml_filename}")
        
        # Save as GML
        gml_filename = filename.replace('.json', '_network.gml')
        nx.write_gml(G, gml_filename)
        print(f"Graph saved as GML: {gml_filename}")
        
        # Print edge type analysis
        print(f"\nDetailed edge type analysis:")
        per_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'per']
        poly_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'poly']
        both_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('edge_type') == 'both']
        
        print(f"- Perfluorinated relationships: {len(per_edges)}")
        print(f"- Polyfluorinated relationships: {len(poly_edges)}")
        print(f"- Both type relationships: {len(both_edges)}")
                
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        print("Make sure the file exists in the current directory.")
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in file '{filename}': {e}")
    except Exception as e:
        print(f"Error: {e}")