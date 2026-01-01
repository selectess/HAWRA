import json
import os
import sys
import argparse

# Add project root to python path
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.append(PROJECT_ROOT)

from bioos.simulations.multiphysics_simulator.simulator import MultiphysicsSimulator

# Configuration de base pour les expériences
base_config = {
    'name': 'pulse_response_experiment',
    'max_time': 200,
    'dt': 0.5,
    'env': {
        'pulse_configs': [
            {'start': 10, 'end': 20, 'intensity': 1.0},
            {'start': 50, 'end': 55, 'intensity': 0.8},
            {'start': 90, 'end': 92, 'intensity': 1.0},
            {'start': 120, 'end': 140, 'intensity': 0.6},
            {'start': 160, 'end': 170, 'intensity': 1.0}
        ]
    },
    'bio': {
        'p700_initial': 0.0,
        'degradation_rate': 0.05,
        'synthesis_rate': 0.2
    },
    'quantum': {
        'threshold': 0.85,
        'decoherence_rate': 0.02
    }
}

def run_experiment(config, output_file):
    """
    Exécute une seule simulation avec une configuration donnée et sauvegarde les résultats.
    """
    print(f"--- Running Experiment: {config.get('name', 'unnamed')} ---")

    output_dir = os.path.dirname(output_file)
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Initialiser et exécuter le simulateur
    simulator = MultiphysicsSimulator(config)
    log = simulator.run()

    # Sauvegarder le log
    log_path = os.path.join(output_dir, "simulation.json")
    with open(log_path, 'w') as f:
        json.dump(log, f, indent=2)
    print(f"Simulation log saved to {log_path}")

    # Générer le graphique final
    simulator.plot_results(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run a pulse-response simulation for HAWRA.')
    parser.add_argument('--degradation_rate', type=float, default=base_config['bio']['degradation_rate'],
                        help='Degradation rate of P700.')
    parser.add_argument('--synthesis_rate', type=float, default=base_config['bio']['synthesis_rate'],
                        help='Synthesis rate of P700.')
    parser.add_argument('--output_file', type=str,
                        required=True,
                        help='Path to save the output plot.')
    args = parser.parse_args()

    # Update config with passed arguments
    base_config['bio']['degradation_rate'] = args.degradation_rate
    base_config['bio']['synthesis_rate'] = args.synthesis_rate

    # Exécuter l'expérience
    run_experiment(base_config, args.output_file)
