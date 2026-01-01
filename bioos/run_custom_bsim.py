import json
import os
import sys

# Ajouter les chemins nécessaires au sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(script_dir, '..')) # Correction du chemin racine
sys_path_to_add = os.path.join(project_root, '03_unified_simulator', 'src')
if sys_path_to_add not in sys.path:
    sys.path.insert(0, sys_path_to_add)

from hawra_simulator.simulator import Simulator # Import corrigé

# Définir les chemins
bsim_path = os.path.join(project_root, 'bioos', 'bio_compiler', 'arbol', 'e2e_validation.bsim.json')

# Initialiser le Simulator avec le chemin du script BSIM
simulator = Simulator(bsim_script=bsim_path)

# Exécuter la simulation
results = simulator.run()

# Sauvegarder les résultats
output_path = os.path.join(project_root, '05_data', 'results', 'validation_results.json')
simulator.save_results(output_path, results)

print("Simulation terminée. Résultats sauvegardés.")