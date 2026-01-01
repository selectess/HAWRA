import os
import subprocess
from flask import Flask, render_template, request, flash
import markdown
import tempfile
from qutip import Bloch, basis, mesolve
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Assurez-vous que le chemin vers votre simulateur est correct
import sys
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.append(project_root)
sys.path.append(os.path.join(project_root, '03_unified_simulator', 'src'))

from unified_simulator import run_simulation_from_bsim

app = Flask(__name__)
app.secret_key = 'super_secret_key' # Replace with a strong secret key in production

@app.route('/')
def index():
    # Pour la démo unifiée, on sert directement le fichier index.html qui se trouve maintenant dans docs/
    # mais Flask attend des templates dans le dossier templates/.
    # Nous allons lire le fichier docs/index.html et le servir.
    docs_index_path = os.path.join(project_root, 'docs', 'index.html')
    with open(docs_index_path, 'r') as f:
        return f.read()

@app.route('/plasmid_3d')
def plasmid_3d():
    # Redirige vers la racine car tout est unifié maintenant
    return index()

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
        env = os.environ.copy()
        env['PYTHONPATH'] = project_root + os.pathsep + os.path.join(project_root, '03_unified_simulator', 'src') + os.pathsep + env.get('PYTHONPATH', '')
        
        compile_command = [
            sys.executable,
            "-m", "arbol.compiler.compiler",
            temp_arbol_path
        ]
        compile_result = subprocess.run(compile_command, capture_output=True, text=True, check=True, env=env, cwd=project_root)
        
        # Run simulation
        simulation_result = run_simulation_from_bsim(bsim_file_path)

        # Generate Bloch sphere plot
        b = Bloch()
        if hasattr(simulation_result, 'states') and simulation_result.states:
            for state in simulation_result.states:
                b.add_states(state)
        
        bloch_output_dir = os.path.join(os.path.dirname(__file__), 'static', 'simulation_results')
        os.makedirs(bloch_output_dir, exist_ok=True)
        filename = f"bloch_{os.path.basename(temp_arbol_path)}.png"
        bloch_image_file = os.path.join(bloch_output_dir, filename)
        b.save(bloch_image_file)
        
        image_file = os.path.join('simulation_results', filename)
        
        # Extract metadata for the UI
        fidelity = 0.95 # Default/Mock if not in result
        if 'fidelity' in simulation_result.data:
            fidelity = simulation_result.data['fidelity']
            
        return render_template('index.html', 
                             arbol_code=arbol_code, 
                             image_file=image_file, 
                             fidelity=fidelity,
                             compilation_output=compile_result.stdout,
                             success=True)

    except subprocess.CalledProcessError as e:
        return render_template('index.html', 
                             arbol_code=arbol_code, 
                             error=f"Compilation failed: {e.stderr}",
                             compilation_output=e.stdout,
                             success=False)
    except Exception as e:
        return render_template('index.html', 
                             arbol_code=arbol_code, 
                             error=f"Simulation error: {str(e)}",
                             success=False)
    finally:
        # Clean up temporary files
        if os.path.exists(temp_arbol_path):
            os.remove(temp_arbol_path)
        if os.path.exists(bsim_file_path):
            os.remove(bsim_file_path)

if __name__ == '__main__':
    app.run(debug=True, port=5005)
