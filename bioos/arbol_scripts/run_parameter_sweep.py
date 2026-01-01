import os
import itertools
from run_pulse_experiment import run_experiment, base_config

def run_parameter_sweep():
    """
    Exécute une série de simulations en faisant varier les paramètres
    pour étudier la sensibilité du modèle.
    """
    print("Starting parameter sweep...")

    # Définir les paramètres à faire varier
    degradation_rates = [0.05, 0.1, 0.2]
    synthesis_rates = [0.2, 0.5, 1.0]

    # Créer le répertoire de base pour les résultats du balayage
    sweep_base_dir = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/arbol_experiments/parameter_sweep"
    if not os.path.exists(sweep_base_dir):
        os.makedirs(sweep_base_dir)

    # Itérer sur toutes les combinaisons de paramètres
    param_combinations = list(itertools.product(degradation_rates, synthesis_rates))
    total_runs = len(param_combinations)
    print(f"Found {total_runs} parameter combinations to run.")

    for i, (deg_rate, syn_rate) in enumerate(param_combinations):
        print(f"\n--- Running simulation {i+1}/{total_runs} ---")
        print(f"Parameters: degradation_rate={deg_rate}, synthesis_rate={syn_rate}")

        # Créer une copie de la configuration de base et la modifier
        config = base_config.copy()
        config['bio'] = config['bio'].copy() # Assurer une copie profonde pour le dict bio
        config['bio']['degradation_rate'] = deg_rate
        config['bio']['synthesis_rate'] = syn_rate

        # Créer un répertoire de sortie unique pour cette simulation
        run_name = f"deg_{deg_rate}_syn_{syn_rate}"
        output_dir = os.path.join(sweep_base_dir, run_name)
        
        # Exécuter l'expérience
        run_experiment(config, output_dir)

    print("\nParameter sweep finished.")

if __name__ == "__main__":
    run_parameter_sweep()
