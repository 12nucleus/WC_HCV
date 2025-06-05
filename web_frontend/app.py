import os
import subprocess
import threading
import queue
from flask import Flask, render_template, request, jsonify, send_from_directory, Response
from pathlib import Path
import datetime
import re
import random # Import random for SVG generation
import math # Import math for SVG generation

app = Flask(__name__)
app.config['SECRET_KEY'] = 'your_secret_key_here' # Replace with a strong secret key in production

# Define base directories relative to the current working directory of the pipeline
# This assumes the Flask app is run from the HCV_pipeline directory or knows its path
PIPELINE_BASE_DIR = Path(os.getcwd())
REPORTS_DIR = PIPELINE_BASE_DIR / "Reports"
CLUSTER_DIRS_PARENT = PIPELINE_BASE_DIR # Cluster_X directories are directly under PIPELINE_BASE_DIR

# Ensure reports directory exists
REPORTS_DIR.mkdir(parents=True, exist_ok=True)

# --- SVG Generation Logic (moved from generate_background_svg.py) ---
VIRUS_PARTICLE_SVG_1 = """
  <!-- Outer envelope -->
  <circle cx="200" cy="200" r="150" fill="#e74c3c" opacity="0.3" stroke="#c0392b" stroke-width="3"/>
  
  <!-- Host-derived lipid envelope -->
  <circle cx="200" cy="200" r="140" fill="none" stroke="#8e44ad" stroke-width="4" stroke-dasharray="8,4"/>
  
  <!-- E1 and E2 envelope glycoproteins (surface spikes) -->
  <g stroke="#27ae60" stroke-width="2" fill="#2ecc71">
    <!-- Multiple spikes around the envelope -->
    <g transform="translate(200,200)">
      <!-- Top spikes -->
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(0)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(30)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(60)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(90)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(120)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(150)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(180)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(210)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(240)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(270)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(300)"/>
      <polygon points="0,-140 -8,-155 8,-155" transform="rotate(330)"/>
    </g>
  </g>
  
  <!-- Nucleocapsid core -->
  <circle cx="200" cy="200" r="100" fill="#3498db" opacity="0.4" stroke="#2980b9" stroke-width="2"/>
  
  <!-- Core protein shell -->
  <circle cx="200" cy="200" r="90" fill="none" stroke="#f39c12" stroke-width="3"/>
  
  <!-- RNA genome (represented as wavy lines inside core) -->
  <g stroke="#e67e22" stroke-width="2" fill="none">
    <path d="M 130 180 Q 150 160, 170 180 T 210 180 T 250 180 T 270 180"/>
    <path d="M 130 200 Q 150 220, 170 200 T 210 200 T 250 200 T 270 200"/>
    <path d="M 130 220 Q 150 200, 170 220 T 210 220 T 250 220 T 270 220"/>
    <path d="M 140 160 Q 160 140, 180 160 T 220 160 T 260 160"/>
    <path d="M 140 240 Q 160 260, 180 240 T 220 240 T 260 240"/>
  </g>
"""

VIRUS_PARTICLE_SVG_2 = """
  <!-- Light gray interior -->
  <ellipse cx="200" cy="200" rx="138" ry="128" fill="#e8e8e8"/>
  
  <!-- Outer envelope with irregular shape -->
  <ellipse cx="200" cy="200" rx="155" ry="145" fill="#e0e0e0" opacity="0.8" stroke="#ee5a52" stroke-width="2"/>
  
  <!-- Host-derived lipid bilayer (double membrane) -->
  <ellipse cx="200" cy="200" rx="145" ry="135" fill="none" stroke="#9c88ff" stroke-width="3"/>
  <ellipse cx="200" cy="200" rx="138" ry="128" fill="none" stroke="#9c88ff" stroke-width="2" opacity="0.6"/>
  
  <!-- E1 and E2 envelope glycoproteins (different arrangement) -->
  <g stroke="#00d2d3" stroke-width="2" fill="#01a3a4">
    <g transform="translate(200,200)">
      <!-- Clustered spikes with varying sizes -->
      <circle cx="0" cy="-135" r="6"/>
      <circle cx="15" cy="-130" r="4"/>
      <circle cx="-12" cy="-132" r="5"/>
      
      <circle cx="95" cy="-95" r="5"/>
      <circle cx="105" cy="-85" r="6"/>
      <circle cx="88" cy="-102" r="4"/>
      
      <circle cx="135" cy="0" r="6"/>
      <circle cx="130" cy="12" r="5"/>
      <circle cx="128" cy="-15" r="4"/>
      
      <circle cx="95" cy="95" r="5"/>
      <circle cx="85" cy="105" r="6"/>
      <circle cx="102" cy="88" r="4"/>
      
      <circle cx="0" cy="135" r="6"/>
      <circle cx="-15" cy="130" r="4"/>
      <circle cx="12" cy="132" r="5"/>
      
      <circle cx="-95" cy="95" r="5"/>
      <circle cx="-105" cy="85" r="6"/>
      <circle cx="-88" cy="102" r="4"/>
      
      <circle cx="-135" cy="0" r="6"/>
      <circle cx="-130" cy="-12" r="5"/>
      <circle cx="-128" cy="15" r="4"/>
      
      <!-- Additional scattered spikes -->
      <circle cx="60" cy="-120" r="3"/>
      <circle cx="-75" cy="-110" r="4"/>
      <circle cx="110" cy="-45" r="3"/>
      <circle cx="120" cy="60" r="4"/>
      <circle cx="45" cy="125" r="3"/>
      <circle cx="-55" cy="115" r="4"/>
      <circle cx="-125" cy="45" r="3"/>
      <circle cx="-115" cy="-65" r="4"/>
    </g>
  </g>
  
  <!-- Nucleocapsid core (hexagonal shape) -->
  <polygon points="200,120 250,150 250,210 200,240 150,210 150,150" fill="#4834d4" opacity="0.5" stroke="#3742fa" stroke-width="2"/>
  
  <!-- Core protein subunits (hexagonal pattern) -->
  <g fill="#ff9ff3" stroke="#f368e0" stroke-width="1">
    <circle cx="200" cy="160" r="6"/>
    <circle cx="180" cy="170" r="6"/>
    <circle cx="220" cy="170" r="6"/>
    <circle cx="180" cy="190" r="6"/>
    <circle cx="220" cy="190" r="6"/>
    <circle cx="200" cy="200" r="6"/>
    <circle cx="165" cy="180" r="5"/>
    <circle cx="235" cy="180" r="5"/>
    <circle cx="200" cy="215" r="5"/>
    <circle cx="185" cy="215" r="5"/>
    <circle cx="215" cy="215" r="5"/>
  </g>
  
  <!-- RNA genome (more complex folded structure) -->
  <g stroke="#ff7675" stroke-width="2.5" fill="none" opacity="0.8">
    <path d="M 165 150 Q 180 135, 200 150 Q 220 165, 235 150 Q 225 170, 210 185"/>
    <path d="M 235 150 Q 250 165, 235 180 Q 220 195, 235 210 Q 220 225, 200 210"/>
    <path d="M 200 210 Q 180 195, 165 210 Q 150 195, 165 180 Q 180 165, 165 150"/>
    <path d="M 185 170 Q 200 155, 215 170 Q 200 185, 185 170"/>
    <path d="M 200 175 Q 215 190, 200 205 Q 185 190, 200 175"/>
  </g>
  
  <!-- Additional membrane details -->
  <g stroke="#a29bfe" stroke-width="1" fill="none" opacity="0.4">
    <ellipse cx="200" cy="200" rx="120" ry="115" stroke-dasharray="3,2"/>
    <ellipse cx="200" cy="200" rx="110" ry="105" stroke-dasharray="2,3"/>
  </g>
"""

VIRUS_PARTICLE_SVG_3 = """
  <!-- Outer envelope (more irregular, slightly flattened) -->
  <ellipse cx="200" cy="200" rx="160" ry="130" fill="#d0d0d0" opacity="0.9" stroke="#8b5a2b" stroke-width="3"/>
  
  <!-- Lipid envelope layers -->
  <ellipse cx="200" cy="200" rx="150" ry="120" fill="none" stroke="#6c5ce7" stroke-width="2" stroke-dasharray="6,3"/>
  <ellipse cx="200" cy="200" rx="140" ry="110" fill="none" stroke="#6c5ce7" stroke-width="2" stroke-dasharray="4,4" opacity="0.7"/>
  
  <!-- E1 and E2 envelope glycoproteins (triangular spikes) -->
  <g stroke="#00b894" stroke-width="1.5" fill="#00cec9">
    <g transform="translate(200,200)">
      <!-- Top region spikes -->
      <polygon points="0,-120 -6,-135 6,-135"/>
      <polygon points="25,-115 19,-130 31,-130"/>
      <polygon points="-20,-118 -26,-133 -14,-133"/>
      <polygon points="45,-100 39,-115 51,-115"/>
      <polygon points="-40,-105 -46,-120 -34,-120"/>
      
      <!-- Right region spikes -->
      <polygon points="140,-10 125,-16 125,-4"/>
      <polygon points="135,15 120,9 120,21"/>
      <polygon points="138,-35 123,-41 123,-29"/>
      <polygon points="132,40 117,34 117,46"/>
      
      <!-- Bottom region spikes -->
      <polygon points="0,120 -6,135 6,135"/>
      <polygon points="-25,115 -31,130 -19,130"/>
      <polygon points="20,118 14,133 26,133"/>
      <polygon points="-45,100 -51,115 -39,115"/>
      <polygon points="40,105 34,120 46,120"/>
      
      <!-- Left region spikes -->
      <polygon points="-140,10 -125,16 -125,4"/>
      <polygon points="-135,-15 -120,-9 -120,-21"/>
      <polygon points="-138,35 -123,41 -123,29"/>
      <polygon points="-132,-40 -117,-34 -117,-46"/>
      
      <!-- Scattered additional spikes -->
      <polygon points="70,-80 64,-95 76,-95"/>
      <polygon points="-65,-85 -71,-100 -59,-100"/>
      <polygon points="85,60 79,75 91,75"/>
      <polygon points="-80,65 -86,80 -74,80"/>
      <polygon points="110,-30 104,-45 116,-45"/>
      <polygon points="-105,35 -111,50 -99,50"/>
    </g>
  </g>
  
  <!-- Nucleocapsid core (circular with segmented appearance) -->
  <circle cx="200" cy="200" r="80" fill="#fd79a8" opacity="0.6" stroke="#e84393" stroke-width="2"/>
  
  <!-- Core protein segments -->
  <g stroke="#2d3436" stroke-width="1" fill="none">
    <circle cx="200" cy="160" r="6"/>
    <circle cx="180" cy="170" r="6"/>
    <circle cx="220" cy="170" r="6"/>
    <circle cx="180" cy="190" r="6"/>
    <circle cx="220" cy="190" r="6"/>
    <circle cx="200" cy="200" r="6"/>
    <circle cx="165" cy="180" r="5"/>
    <circle cx="235" cy="180" r="5"/>
    <circle cx="200" cy="215" r="5"/>
    <circle cx="185" cy="215" r="5"/>
    <circle cx="215" cy="215" r="5"/>
  </g>
  
  <!-- Core protein units (pentagonal shapes) -->
  <g fill="#fdcb6e" stroke="#e17055" stroke-width="1">
    <polygon points="200,140 210,150 205,165 195,165 190,150"/>
    <polygon points="240,160 250,170 245,185 235,185 230,170"/>
    <polygon points="260,200 270,210 265,225 255,225 250,210"/>
    <polygon points="240,240 250,250 245,265 235,265 230,250"/>
    <polygon points="200,260 210,270 205,285 195,285 190,270"/>
    <polygon points="160,240 170,250 165,265 155,265 150,250"/>
    <polygon points="140,200 150,210 145,225 135,225 130,210"/>
    <polygon points="160,160 170,170 165,185 155,185 150,170"/>
  </g>
  
  <!-- RNA genome (spiral pattern) -->
  <g stroke="#e17055" stroke-width="3" fill="none" opacity="0.9">
    <path d="M 200 160 
             Q 220 170, 230 190 
             Q 240 210, 230 230 
             Q 220 250, 200 240 
             Q 180 230, 170 210 
             Q 160 190, 170 170 
             Q 180 150, 200 160"/>
    <path d="M 200 180 
             Q 210 185, 215 195 
             Q 220 205, 215 215 
             Q 210 225, 200 220 
             Q 190 215, 185 205 
             Q 180 195, 185 185 
             Q 190 175, 200 180"/>
  </g>
  
  <!-- Additional membrane texture -->
  <g stroke="#74b9ff" stroke-width="0.8" fill="none" opacity="0.3">
    <ellipse cx="200" cy="200" rx="125" ry="95" stroke-dasharray="2,1"/>
    <ellipse cx="200" cy="200" rx="115" ry="85" stroke-dasharray="1,2"/>
    <ellipse cx="200" cy="200" rx="105" ry="75" stroke-dasharray="3,1"/>
  </g>
"""

VIRUS_PARTICLE_SVGS = [VIRUS_PARTICLE_SVG_1, VIRUS_PARTICLE_SVG_2, VIRUS_PARTICLE_SVG_3]

def generate_random_virus_background(
    canvas_width=800, # Smaller canvas
    canvas_height=450, # Smaller canvas
    num_particles=25, # Increased number of particles for more density
    min_scale=0.1, # Reduced by half
    max_scale=0.2, # Reduced by half
    opacity=0.85 # Set opacity to 85%
):
    svg_elements = []
    for i in range(num_particles):
        x = random.uniform(0, canvas_width)
        y = random.uniform(0, canvas_height)
        scale = random.uniform(min_scale, max_scale)
        rotation = random.uniform(0, 360)

        particle_content = random.choice(VIRUS_PARTICLE_SVGS) # Randomly choose between the two SVGs

        svg_elements.append(f"""
        <g transform="translate({x-200*scale} {y-200*scale}) scale({scale}) rotate({rotation} 200 200)" opacity="{opacity}"
           >
            {particle_content}
        </g>
        """)

    final_svg = f"""
<svg width="{canvas_width}" height="{canvas_height}" viewBox="0 0 {canvas_width} {canvas_height}" xmlns="http://www.w3.org/2000/svg">
    <defs>
        <!-- Define the virus particle as a symbol if needed, but direct embedding is fine for a few instances -->
    </defs>
    {''.join(svg_elements)}
</svg>
"""
    return final_svg.replace('\n', '').replace('\t', '').strip()

def generate_random_virus_background(
    canvas_width=800, # Smaller canvas
    canvas_height=450, # Smaller canvas
    num_particles=25, # Increased number of particles for more density
    min_scale=0.1, # Reduced by half
    max_scale=0.2, # Reduced by half
    opacity=0.85 # Set opacity to 85%
):
    svg_elements = []
    for i in range(num_particles):
        x = random.uniform(0, canvas_width)
        y = random.uniform(0, canvas_height)
        scale = random.uniform(min_scale, max_scale)
        rotation = random.uniform(0, 360)

        # particle_content is now chosen randomly from VIRUS_PARTICLE_SVGS
        selected_particle_content = random.choice(VIRUS_PARTICLE_SVGS)

        svg_elements.append(f"""
        <g transform="translate({x-200*scale} {y-200*scale}) scale({scale}) rotate({rotation} 200 200)" opacity="{opacity}"
           >
            {selected_particle_content}
        </g>
        """)

    final_svg = f"""
<svg width="{canvas_width}" height="{canvas_height}" viewBox="0 0 {canvas_width} {canvas_height}" xmlns="http://www.w3.org/2000/svg">
    <defs>
        <!-- Define the virus particle as a symbol if needed, but direct embedding is fine for a few instances -->
    </defs>
    {''.join(svg_elements)}
</svg>
"""
    return final_svg.replace('\n', '').replace('\t', '').strip()

@app.route('/')
def index():
    generated_svg = generate_random_virus_background()
    return render_template('index.html', background_svg=generated_svg)

@app.route('/run_pipeline', methods=['POST'])
def run_pipeline():
    data = request.json
    reads_dir_input = data.get('reads_dir', 'test4/')
    keep_tmp = data.get('keep_tmp', False)
    keep_unmerged = data.get('keep_unmerged', False)
    r1_pattern = data.get('r1_pattern', '*_R1_001.fastq.gz')
    kmer_overlap_threshold = data.get('kmer_overlap_threshold', 0.005)
    snp_distance_threshold = data.get('snp_distance_threshold', 9)

    # Construct the command for HCV_Master_script_V0_1.py
    # The base_dir argument for the pipeline script is the current working directory
    # where the pipeline scripts (HCV_Master_script_V0_1.py) reside.
    command = [
        'python',
        str(PIPELINE_BASE_DIR / 'HCV_Master_script_V0_1.py'),
        str(PIPELINE_BASE_DIR / reads_dir_input), # reads_dir argument
        str(PIPELINE_BASE_DIR) # base_dir argument
    ]

    if keep_tmp:
        command.append('--keep_tmp')
    if keep_unmerged:
        command.append('--keep_unmerged')
    
    command.extend([
        '--r1_pattern', r1_pattern,
        '--kmer_overlap_threshold', str(kmer_overlap_threshold),
        '--snp_distance_threshold', str(snp_distance_threshold)
    ])

    # Hardcode other script paths and executable paths as they are relative to base_dir
    # and less likely to change per run.
    command.extend([
        '--combine_script', 'HCV_combine_reads_V0_1.py',
        '--distancer_script', 'HCV_kmers_distancer_V0_1.py',
        '--transmission_script', 'HCV_transmission_test_V0_2.py', # Updated from V0_1 to V0_2 based on file list
        '--snp_dists_path', 'snp-dists',
        '--mafft_path', 'mafft',
        '--rscript_path', 'Rscript'
    ])

    def generate():
        process = None
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, # Combine stdout and stderr
                text=True,
                bufsize=1, # Line-buffered
                universal_newlines=True,
                cwd=PIPELINE_BASE_DIR # Run the command from the pipeline's base directory
            )

            for line in iter(process.stdout.readline, ''):
                yield f"data: {line}\n\n"
            
            process.stdout.close()
            return_code = process.wait()
            yield f"data: --- Pipeline Finished with exit code {return_code} ---\n\n"

        except FileNotFoundError:
            yield "data: Error: Python executable or pipeline script not found. Ensure Python is in PATH and script exists.\n\n"
        except Exception as e:
            yield f"data: An unexpected error occurred: {e}\n\n"
        finally:
            if process and process.poll() is None: # If process is still running
                process.terminate()
                process.wait()
                yield "data: --- Pipeline process terminated unexpectedly ---\n\n"

    return Response(generate(), mimetype='text/event-stream')

@app.route('/list_files')
def list_files():
    report_files = []
    html_files = []

    # List .txt files in Reports/
    for f in REPORTS_DIR.glob('*.txt'):
        report_files.append(str(f.relative_to(PIPELINE_BASE_DIR)))
    for f in REPORTS_DIR.glob('*.out'): # Also include .out files as they are reports
        report_files.append(str(f.relative_to(PIPELINE_BASE_DIR)))

    # List .html files in Reports/ and Cluster_X/ directories
    # Search directly under PIPELINE_BASE_DIR for Cluster_X/
    for cluster_dir in CLUSTER_DIRS_PARENT.iterdir():
        if cluster_dir.is_dir() and re.match(r'Cluster_\d+', cluster_dir.name):
            for f in cluster_dir.glob('*.html'):
                html_files.append(str(f.relative_to(PIPELINE_BASE_DIR)))
    
    # Also check Reports/ for HTML files
    for f in REPORTS_DIR.glob('*.html'):
        html_files.append(str(f.relative_to(PIPELINE_BASE_DIR)))

    return jsonify(text_files=sorted(report_files), html_files=sorted(html_files))

@app.route('/view_file/<path:filename>')
def view_file(filename):
    file_path = PIPELINE_BASE_DIR / filename
    
    # Security check: ensure file is within allowed directories
    allowed_base_paths = [REPORTS_DIR]
    for cluster_dir_name in os.listdir(PIPELINE_BASE_DIR):
        if re.match(r'Cluster_\d+', cluster_dir_name):
            allowed_base_paths.append(PIPELINE_BASE_DIR / cluster_dir_name)

    is_safe_path = False
    for base_path in allowed_base_paths:
        try:
            file_path.relative_to(base_path) # This will raise ValueError if not relative
            is_safe_path = True
            break
        except ValueError:
            pass

    if not is_safe_path:
        return "Access Denied", 403

    if not file_path.is_file():
        return "File not found", 404

    if file_path.suffix == '.html':
        # For HTML files, serve them directly
        return send_from_directory(file_path.parent, file_path.name)
    else:
        # For text files, read content and return
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                content = f.read()
            return Response(content, mimetype='text/plain')
        except Exception as e:
            return f"Error reading file: {e}", 500

@app.route('/download_file/<path:filename>')
def download_file(filename):
    file_path = PIPELINE_BASE_DIR / filename

    # Security check: ensure file is within allowed directories
    allowed_base_paths = [REPORTS_DIR]
    for cluster_dir_name in os.listdir(PIPELINE_BASE_DIR):
        if re.match(r'Cluster_\d+', cluster_dir_name):
            allowed_base_paths.append(PIPELINE_BASE_DIR / cluster_dir_name)

    is_safe_path = False
    for base_path in allowed_base_paths:
        try:
            file_path.relative_to(base_path) # This will raise ValueError if not relative
            is_safe_path = True
            break
        except ValueError:
            pass

    if not is_safe_path:
        return "Access Denied", 403

    if not file_path.is_file():
        return "File not found", 404

    return send_from_directory(file_path.parent, file_path.name, as_attachment=True)

if __name__ == '__main__':
    # In a production environment, use a more robust WSGI server like Gunicorn or uWSGI
    app.run(debug=True, host='0.0.0.0', port=5000)
