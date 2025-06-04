#!/usr/bin/env python3
import networkx as nx
import numpy as np
from collections import defaultdict
from pyvis.network import Network
from PIL import Image
from scipy.sparse.csgraph import minimum_spanning_tree
import json
import matplotlib
import matplotlib.pyplot as plt
import sys
import colorsys
import hashlib
import io
import base64
import re # Import regex module

import argparse # Added for CLI argument parsing

matplotlib.use('Agg')
target_sample=""

def read_distance_matrix(filename):
    samples = []
    dist_matrix = []
    
    with open(filename, 'r') as f:
        # Read header line to get sample names
        samples = f.readline().strip().split('\t')
        
        # Read distance matrix
        for line in f:
            row = [float(x) for x in line.strip().split('\t')[1:]]  # Skip sample name
            dist_matrix.append(row)
    
    return samples, dist_matrix

def build_eburst_mst(samples, dist_matrix):
    # Create lookup for sample indices
    sample_to_idx = {sample: idx for idx, sample in enumerate(samples)}
    
    # Track which samples have been assigned to clusters
    assigned_samples = set()
    clusters = {}
    cluster_id = 0
    
    # Group samples into clusters (only 0 SNP differences)
    for i in range(len(samples)):
        if samples[i] in assigned_samples:
            continue
            
        current_cluster = set([samples[i]])
        assigned_samples.add(samples[i])
        
        for j in range(len(samples)):
            if i != j and samples[j] not in assigned_samples:
                # Check if sample has 0 SNP difference with ALL current cluster members
                is_identical = True
                for member in current_cluster:
                    member_idx = sample_to_idx[member]
                    if dist_matrix[sample_to_idx[samples[j]]][member_idx] != 0:
                        is_identical = False
                        break
                
                if is_identical:
                    current_cluster.add(samples[j])
                    assigned_samples.add(samples[j])
        
        if current_cluster:
            clusters[cluster_id] = list(current_cluster)
            cluster_id += 1
    
    # Create MST between clusters
    G = nx.Graph()
    for i in range(len(clusters)):
        G.add_node(i)
    
    # Calculate average distances between clusters
    for i in range(len(clusters)):
        for j in range(i + 1, len(clusters)):
            distances = []
            for sample1 in clusters[i]:
                for sample2 in clusters[j]:
                    idx1 = sample_to_idx[sample1]
                    idx2 = sample_to_idx[sample2]
                    distances.append(dist_matrix[idx1][idx2])
            if distances: # Avoid division by zero if a cluster is somehow empty (shouldn't happen with current logic)
                avg_dist = sum(distances) / len(distances)
                G.add_edge(i, j, weight=avg_dist)
    
    # Handle case where G might be empty or disconnected
    if not G.nodes():
        return nx.Graph(), {} # Return empty graph and clusters if no nodes were added
        
    # Compute MST - handle potential errors if graph is not connected
    try:
        mst = nx.minimum_spanning_tree(G)
    except nx.NetworkXNoCycle:
        # If no cycles, the graph itself is a forest (or disconnected nodes)
        # The MST is the graph itself in this case
        mst = G
    except Exception as e:
        print(f"Error computing MST: {e}", file=sys.stderr)
        # Depending on desired behavior, might return empty graph or raise error
        return nx.Graph(), clusters # Return empty graph but keep clusters


    return mst, clusters

def save_cluster_info(mst, clusters, filename="clusters.txt"):
    with open(filename, 'w') as f:
        for cluster_id, samples in clusters.items():
            # Write cluster header and samples
            f.write(f"Node{cluster_id} ({len(samples)} samples):\n")
            for sample in samples:
                f.write(f"  {sample}\n")
            
            # Get neighboring clusters and distances in the MST
            neighbors = mst.adj.get(cluster_id, {}) # Use .adj.get for safety if node not in MST
            if neighbors:
                f.write("\nNeighboring clusters (in MST):\n")
                for neighbor, edge_data in neighbors.items():
                    distance = edge_data['weight']
                    f.write(f"  Node {neighbor}: {distance:.0f} SNPs\n")
            f.write("\n")

def generate_gradient_color(index, total):
    """Generate a color gradient from light gray -> dark blue -> dark red -> light orange."""
    if total == 0: # Handle case with no unique sample types
        return '#CCCCCC' # Default gray
    if total == 1:
        return '#0000FF'  # Dark blue for single sample type
    
    # Define color stops
    colors = [
        (220, 220, 220),  # Light gray
        (0, 0, 150),      # Dark blue
        (150, 0, 0),      # Dark red
        (255, 200, 100)   # Light orange
    ]
    
    # Calculate position in gradient (0 to 3 for 4 color stops)
    position = (index / (total - 1)) * (len(colors) - 1)
    
    # Find the two colors to interpolate between
    color_index = int(position)
    color_blend = position - color_index
    
    if color_index >= len(colors) - 1:
        return '#{:02x}{:02x}{:02x}'.format(*colors[-1])
    
    # Interpolate between the two colors
    c1 = colors[color_index]
    c2 = colors[color_index + 1]
    
    r = int(c1[0] * (1 - color_blend) + c2[0] * color_blend)
    g = int(c1[1] * (1 - color_blend) + c2[1] * color_blend)
    b = int(c1[2] * (1 - color_blend) + c2[2] * color_blend)
    
    return '#{:02x}{:02x}{:02x}'.format(r, g, b)

def create_pie_chart_node(sample_counts, sample_colors, node_id, node_size=50):
    """Create a pie chart for a node with multiple sample types and return as base64 encoded image."""
    # Create figure with transparent background
    fig, ax = plt.subplots(figsize=(1, 1), dpi=200)
    fig.patch.set_alpha(0.0)
    ax.set_aspect('equal')
    
    # Get sample types and their counts
    labels = list(sample_counts.keys())
    sizes = list(sample_counts.values())
    
    # Ensure colors match the labels from sample_counts
    colors = [sample_colors.get(label, '#CCCCCC') for label in labels] # Use .get with default for safety
    
    # Create pie chart
    ax.pie(sizes, colors=colors, startangle=90, wedgeprops={'edgecolor': 'black', 'linewidth': 0.5})
    ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle
    
    # Convert plot to image
    buf = io.BytesIO()
    fig.savefig(buf, format='png', transparent=True, bbox_inches='tight', pad_inches=0)
    plt.close(fig)
    buf.seek(0)
    
    # Convert to base64
    img_str = base64.b64encode(buf.read()).decode('ascii')
    return f"data:image/png;base64,{img_str}"


def parse_sample_name(sample_string):

    parts = sample_string.split('_')
    if len(parts) <= 2:
        return sample_string
    return '_'.join(parts[1:-1])

def plot_interactive_mst(mst, clusters, output_file, my_min_degree,min_snp):
    net = Network(height="750px", width="100%", bgcolor="#ffffff", font_color="black")
    
    # Track unique base sample names for coloring
    base_sample_names_for_coloring = set()
    # Also track all unique base sample names to ensure representation
    all_unique_base_samples = set()

    # First pass - collect all unique base sample names for coloring and for ensuring representation
    for node_id, samples_in_cluster in clusters.items(): # Iterate through clusters directly
        print(f"xxxxx Processing cluster {node_id} with {samples_in_cluster} samples")
        for sample in samples_in_cluster:
            base_name = parse_sample_name(sample)
            base_sample_names_for_coloring.add(base_name)
            all_unique_base_samples.add(base_name) # Collect all unique base samples

    # Create color mapping based on unique base sample names
    base_sample_names_list = sorted(list(base_sample_names_for_coloring))
    sample_colors = {
        base_name: generate_gradient_color(i, len(base_sample_names_list))
        for i, base_name in enumerate(base_sample_names_list)
    }
    
    # Calculate degrees in the MST
    degrees = dict(mst.degree())

    # --- New Logic for ensuring sample representation ---
    guaranteed_represented_samples = set()
    nodes_to_force_include = set() # Nodes that must be included to represent a sample

    # First pass: Identify nodes that must be included to represent unrepresented samples
    # Sort nodes to ensure consistent behavior if multiple nodes could represent the same sample
    # (e.g., by node ID or number of samples in cluster)
    sorted_nodes = sorted(mst.nodes(), key=lambda n: len(clusters.get(n, []))) # Sort by number of samples, smaller first

    for node in sorted_nodes:
        samples_in_node = clusters.get(node, [])
        node_base_samples = {parse_sample_name(s) for s in samples_in_node}
        
        # Check if this node represents any base sample not yet guaranteed
        newly_represented_samples = node_base_samples - guaranteed_represented_samples
        
        if newly_represented_samples:
            nodes_to_force_include.add(node)
            guaranteed_represented_samples.update(newly_represented_samples)
            
            # Optimization: If all samples are now guaranteed, we can stop this pass
            if guaranteed_represented_samples == all_unique_base_samples:
                break
    # --- End New Logic ---

    # Add nodes - Filter nodes with degree <= 1 OR if not forcefully included
    visible_nodes = set()
    for node in mst.nodes():
        # --- Filtering: Only add nodes with degree > 1 in the MST OR if forcefully included ---
        if node in nodes_to_force_include:
            # This node is forcefully included to ensure sample representation
            pass # Do not apply degree filter
        elif degrees.get(node, 0) <= my_min_degree:
            # Original filtering: skip if degree is too low
            continue

        samples = clusters.get(node, [])
        base_name_counts = {} # Count occurrences of each base sample name
        
        for sample in samples:
            base_name = parse_sample_name(sample) # Use the revised parse_sample_name
            
            # Count base sample names for pie chart
            base_name_counts[base_name] = base_name_counts.get(base_name, 0) + 1
        
        node_size = 10 * len(samples) # Node size based on total samples in cluster
        
        # Determine node color or create pie chart based on base sample names
        if len(base_name_counts) <= 1 or sum(base_name_counts.values()) == 0:
            # Single base sample name type or no valid samples - use solid color
            if base_name_counts:
                dominant_base_name = list(base_name_counts.keys())[0]
                node_color = sample_colors.get(dominant_base_name, '#CCCCCC')
            else:
                node_color = '#CCCCCC'
                
            net.add_node(
                node,
                label=f"Node {node}\n({len(samples)} samples)",
                title="\n".join(samples), # Keep original sample names in title
                size=node_size,
                color=node_color
            )
        else:
            # Multiple base sample names - use pie chart
            pie_chart_data = create_pie_chart_node(base_name_counts, sample_colors, node, node_size)
            
            net.add_node(
                node,
                label=f"Node {node}\n({len(samples)} samples)",
                title="\n".join(samples), # Keep original sample names in title
                size=node_size,
                shape="image",
                image=pie_chart_data,
                sampleCount=len(samples)
            )
        visible_nodes.add(node)
    
    # Add edges - Filter edges connected to nodes with degree <= 1 AND edges with weight > 9
    for (u, v, d) in mst.edges(data=True):
        # Only add edges if both connected nodes are visible AND weight <= 9
        if u in visible_nodes and v in visible_nodes and d['weight'] <= 9:
            net.add_edge(
                u, v,
                label=f"{d['weight']:.0f}",
                title=f"{d['weight']:.0f} SNPs",
                length=200
            )
    
    # Updated physics configuration
    physics_config = {
        "physics": {
            "enabled": True,
            "barnesHut": {
                "gravitationalConstant": -2000,
                "centralGravity": 0.3,
                "springLength": 200,
                "springConstant": 0.04
            },
            "minVelocity": 0.75,
            "solver": "barnesHut"
        }
    }
    
    net.set_options(json.dumps(physics_config))
    
    # Create legend HTML for sample types (now base sample names)
    legend_html = """
    <div style="position: fixed; top: 10px; right: 10px; background: white;
                padding: 10px; border: 1px solid black; border-radius: 5px; z-index: 1000;">
        <h3 style="margin-top: 0;">Base Sample Names</h3>
        <ul style="list-style-type: none; padding: 0; margin: 0;">
    """
    
    # Add entries for each base sample name
    for base_name, color in sorted(sample_colors.items()):
        legend_html += f"""
            <li style="margin: 5px; display: flex; align-items: center;">
                <span style="display: inline-block; width: 20px; height: 20px;
                      background-color: {color}; margin-right: 5px; border: 1px solid black;"></span>
                {base_name}
            </li>
        """
    
    legend_html += f"""
        </ul>
        <div style="margin-top: 10px; font-style: italic;">
            Note: Only showing edges with ≤ {min_snp} SNPs distance<br>
            Only showing nodes with > {my_min_degree} connection in MST
        </div>
    </div>
    """
    
    # Create controls HTML with improved export functionality
    controls_html = """
    <div style="position: fixed; bottom: 10px; left: 10px; background: white;
                padding: 10px; border: 1px solid black; border-radius: 5px; z-index: 1000;">
        <button id="playPauseBtn">Pause</button>
        <button id="savePngBtn">Save to PNG</button>
        <button id="savePdfBtn">Save to PDF</button>
        <div style="margin-top: 10px;">
            <label>Gravity: </label>
            <input type="range" id="gravitySlider" min="-15000" max="0" value="-2000">
            <span id="gravityValue">-2000</span>
        </div>
        <div>
            <label>Spring Length: </label>
            <input type="range" id="springLengthSlider" min="50" max="500" value="200">
            <span id="springLengthValue">200</span>
        </div>
        <div>
            <label>Spring Constant: </label>
            <input type="range" id="springConstantSlider" min="0" max="0.5" value="0.04" step="0.01">
            <span id="springConstantValue">0.04</span>
        </div>
    </div>
    
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js"></script>
    <script>
    function initializeControls() {
        if (!window.network) { 
            console.error('Network not initialized'); 
            return; 
        }
        
        document.getElementById('playPauseBtn').addEventListener('click', function() {
            var enabled = window.network.physics.options.enabled;
            window.network.setOptions({physics:{enabled:!enabled}});
            this.innerText = enabled ? 'Play' : 'Pause';
        });
        
        document.getElementById('savePngBtn').addEventListener('click', function() {
            saveCanvasAsImage('png');
        });
        
        document.getElementById('savePdfBtn').addEventListener('click', function() {
            saveCanvasAsImage('pdf');
        });
        
        ['gravity','springLength','springConstant'].forEach(function(param) {
            var slider = document.getElementById(param + 'Slider');
            slider.addEventListener('input', function() { 
                updatePhysics(param, this.value); 
            });
        });
        
        console.log("Controls initialized successfully");
    }
    
    function saveCanvasAsImage(format) {
        var canvas = document.querySelector('#mynetwork canvas');
        if (!canvas) {
            console.error('Canvas not found');
            return;
        }
        
        // Create a new canvas with white background
        var newCanvas = document.createElement('canvas');
        
        // Set high resolution (10x10 inches at 300 DPI)
        var dpi = 300;
        newCanvas.width = 5 * dpi;
        newCanvas.height = 5 * dpi;
        
        var context = newCanvas.getContext('2d');
        
        // Fill with white background
        context.fillStyle = 'white';
        context.fillRect(0, 0, newCanvas.width, newCanvas.height);
        
        // Scale the drawing to fit the higher resolution canvas
        var scale = Math.min(
            newCanvas.width / canvas.width,
            newCanvas.height / canvas.height
        );
        
        // Center the drawing
        var offsetX = (newCanvas.width - canvas.width * scale);
        var offsetY = (newCanvas.height - canvas.height * scale);
        
        // Draw the original canvas content on top with scaling
        context.save();
        context.translate(offsetX, offsetY);
        context.scale(scale, scale);
        context.drawImage(canvas, 0, 0);
        context.restore();
        
        // Default filename
        var defaultFilename = 'mst_visualization';
        var filename = prompt('Save visualization as:', defaultFilename);
        if (!filename) return;
        
        if (format === 'png') {
            // For PNG format
            if (!filename.toLowerCase().endsWith('.png')) {
                filename += '.png';
            }
            
            newCanvas.toBlob(function(blob) {
                var a = document.createElement('a');
                a.href = URL.createObjectURL(blob);
                a.download = filename;
                a.style.display = 'none';
                document.body.appendChild(a);
                a.click();
                setTimeout(function() {
                    document.body.removeChild(a);
                    URL.revokeObjectURL(a.href);
                }, 100);
            }, 'image/png');
        } else if (format === 'pdf') {
            // For PDF format
            if (!filename.toLowerCase().endsWith('.pdf')) {
                filename += '.pdf';
            }
            
            try {
                // Create PDF with jsPDF
                var pdf = new jspdf.jsPDF({
                    orientation: 'landscape',
                    unit: 'in',
                    format: [10, 10]
                });
                
                // Add the image to the PDF
                var imgData = newCanvas.toDataURL('image/png');
                pdf.addImage(imgData, 'PNG', 0, 0, 10, 10);
                
                // Save the PDF
                pdf.save(filename);
            } catch (e) {
                console.error('Error creating PDF:', e);
                alert('Error creating PDF. Check console for details.');
            }
        }
    }
    
    
    
    function updatePhysics(param, value) {
        console.log("Inside updatePhysics function for param:", param, "value:", value);
        if (!window.network) {
             console.error("updatePhysics Error: window.network is not defined!");
             return;
        }
        
        // Create minimal physics options with just the parameter we want to change
        var physicsOptions = {
             physics: {
                  barnesHut: {}
             }
        };

        if (param === 'gravity') {
            physicsOptions.physics.barnesHut.gravitationalConstant = parseInt(value);
            document.getElementById('gravityValue').innerText = value;
        } else if (param === 'springLength') {
            physicsOptions.physics.barnesHut.springLength = parseInt(value);
            document.getElementById('springLengthValue').innerText = value;
        } else if (param === 'springConstant') {
            physicsOptions.physics.barnesHut.springConstant = parseFloat(value);
            document.getElementById('springConstantValue').innerText = value;
        } else {
             console.warn("updatePhysics: Unknown parameter:", param);
             return;
        }
        
        console.log("Setting physics options:", JSON.stringify(physicsOptions, null, 2));
        try {
             window.network.setOptions(physicsOptions);
             // Ensure physics stays enabled and button text is correct
             if (!window.network.physics.options.enabled) {
                  window.network.setOptions({ physics: { enabled: true } });
             }
             document.getElementById('playPauseBtn').innerText = 'Pause';
             console.log("Physics options updated successfully.");
        } catch (e) {
             console.error("Error setting physics options:", e);
        }
    }
    
    // Enhanced initialization with debug logging
    function initVisualization() {
        console.log("Initializing visualization controls...");
        
        // Get reference to network object
        window.network = network;
        
        if (!window.network) {
            console.error("Network object not found!");
            return;
        }
        
        // Initialize all controls
        try {
            initializeControls();
            
            // Add toggle button handler (if you add one later)
            // document.getElementById('toggleSingleBtn').addEventListener('click', function() {
            //     console.log("Toggle button clicked");
            //     toggleSingleSamples.call(this);
            // });
            
            // Add physics control handlers
            ['gravity','springLength','springConstant'].forEach(function(param) {
                var slider = document.getElementById(param + 'Slider');
                if (slider) {
                    slider.addEventListener('input', function() {
                        console.log(param + " slider changed to " + this.value);
                        updatePhysics(param, this.value);
                    });
                } else {
                    console.error(param + " slider not found!");
                }
            });
            
            console.log("All controls initialized successfully");
        } catch (e) {
            console.error("Error initializing controls:", e);
        }
    }
    
    // Initialize when both DOM and network are ready
    document.addEventListener('DOMContentLoaded', function() {
        console.log("DOM content loaded");
        
        if (typeof network !== 'undefined') {
            initVisualization();
        } else {
            console.log("Waiting for network initialization...");
            var checkInterval = setInterval(function() {
                if (typeof network !== 'undefined') {
                    clearInterval(checkInterval);
                    initVisualization();
                }
            }, 500);
        }
    });
    </script>
    """
    
    # Save graph
    net.save_graph(output_file)
    
    # Read the file and inject both legend and controls
    with open(output_file, 'r') as f:
        content = f.read()
    
    # Insert controls and legend before the closing body tag
    content = content.replace('</body>', f'{controls_html}{legend_html}</body>')
    
    # Update the drawGraph function to export network to window object
    content = content.replace('return network;', 'window.network = network;\nreturn network;')
    
    # Add jsPDF library for PDF export (already present, but ensure it's not duplicated)
    # content = content.replace('<head>', '<head>\n<script src="https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js"></script>')
    
    with open(output_file, 'w') as f:
        f.write(content)

def main():
    parser = argparse.ArgumentParser(description="Builds an interactive Minimum Spanning Tree (MST) visualization from a SNP distance matrix.")
    parser.add_argument("input_file", type=str, help="Path to the input SNP distance matrix file (TSV format).")
    parser.add_argument("output_file", type=str, help="Path to the output HTML file for the interactive MST visualization.")
    parser.add_argument("--min_degree", type=int, default=1,help="Minimum degree for filtering nodes in the visualization. Nodes with fewer connections than this value will be hidden unless they are required to represent a unique sample type. Default is 1.")
    parser.add_argument("--min_snp", type=int, default=9, help="Minimum SNP distance for edges to be included in the visualization. Edges with a weight greater than this value will be excluded. Default is 9.")
    args = parser.parse_args()

    global my_min_degree
    global my_min_snp
    my_min_degree = args.min_degree # Assign from parsed arguments
    my_min_snp = args.min_snp # Assign from parsed arguments
    samples, dist_matrix = read_distance_matrix(args.input_file)
    mst, clusters = build_eburst_mst(samples, dist_matrix)
    save_cluster_info(mst, clusters)
    plot_interactive_mst(mst, clusters, args.output_file, my_min_degree, my_min_snp)
    print(f"MST visualization generated successfully. Open {args.output_file} in your web browser.")

if __name__ == "__main__":
    main()

# Install required packages:
# pip install pyvis networkx scipy matplotlib pillow
