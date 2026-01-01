import os
import subprocess
from flask import Flask, render_template, request, flash
import markdown
import tempfile
from qutip import Bloch, basis, mesolve
import numpy as np
import matplotlib.pyplot as plt

# Assurez-vous que le chemin vers votre simulateur est correct
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '03_unified_simulator', 'src')))
from unified_simulator import run_simulation_from_bsim

app = Flask(__name__)
app.secret_key = 'super_secret_key' # Replace with a strong secret key in production

@app.route('/')
def index():
    presentation_path = os.path.join(os.path.dirname(__file__), '../06_presentation/presentation.md')
    with open(presentation_path, 'r') as f:
        md_content = f.read()
    html_content = markdown.markdown(md_content)
    return render_template('index.html', content=html_content)

@app.route('/plasmid_3d')
def plasmid_3d():
    return render_template('plasmid_3d.html')

@app.route('/compile_and_simulate', methods=['POST'])
def compile_and_simulate():
    arbol_code = request.form['arbol_code']
    
    # Create a temporary Arbol file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.arbol') as temp_arbol_file:
        temp_arbol_file.write(arbol_code)
        temp_arbol_path = temp_arbol_file.name

    bsim_file_path = temp_arbol_path.replace('.arbol', '.bsim.json')
    
    try:
        # Compile Arbol code
        compile_command = [
            'python',
            os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'arbol', 'compiler', 'compiler.py')),
            temp_arbol_path
        ]
        compile_result = subprocess.run(compile_command, capture_output=True, text=True, check=True)
        flash(f"Compilation successful: {compile_result.stdout}", 'success')

        # Run simulation
        simulation_results = run_simulation_from_bsim(bsim_file_path)

        # Generate Bloch sphere plot
        b = Bloch()
        # Assuming simulation_results contains states that can be plotted on Bloch sphere
        # This part might need adjustment based on the actual output of run_simulation_from_bsim
        # For example, if it returns a list of Qobj states:
        for state in simulation_results.states:
            b.add_states(state)
        
        bloch_output_dir = os.path.join(os.path.dirname(__file__), 'static', 'simulation_results')
        os.makedirs(bloch_output_dir, exist_ok=True)
        bloch_image_file = os.path.join(bloch_output_dir, 'bloch_sphere.png')
        b.save(bloch_image_file)
        
        flash("Simulation successful and Bloch sphere generated.", 'info')
        image_file = os.path.join('simulation_results', 'bloch_sphere.png')

    except subprocess.CalledProcessError as e:
        flash(f"Compilation failed: {e.stderr}", 'error')
        image_file = None
    except Exception as e:
        flash(f"An unexpected error occurred during simulation or Bloch sphere generation: {e}", 'error')
        image_file = None
    finally:
        # Clean up temporary files
        os.remove(temp_arbol_path)
        if os.path.exists(bsim_file_path):
            os.remove(bsim_file_path)

    presentation_path = os.path.join(os.path.dirname(__file__), '../06_presentation/presentation.md')
    with open(presentation_path, 'r') as f:
        md_content = f.read()
    html_content = markdown.markdown(md_content)

    return render_template('index.html', content=html_content, arbol_code=arbol_code, image_file=image_file)

if __name__ == '__main__':
    app.run(debug=True, port=8080)
