
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, sigmax, mesolve

# Importer les modèles précédents
from cry2_model import simulate_cry2_expression
from metabolic_model import simulate_metabolism

def run_multiphysics_simulation(em_frequency, sim_duration=10):
    """
    Exécute une simulation multiphysique intégrée.

    Args:
        em_frequency (float): Fréquence du signal EM en Hz.
        sim_duration (float): Durée de la simulation en secondes.

    Returns:
        dict: Un dictionnaire contenant les résultats de la simulation.
    """
    # 1. Simulation biochimique (CRY2)
    _, cry2_expression_over_time = simulate_cry2_expression(em_frequency, sim_duration)
    final_cry2_level = cry2_expression_over_time[-1]

    # 2. Simulation quantique (Bio-qubit P700)
    # L'angle de rotation de la porte quantique dépend du niveau d'expression de CRY2
    # C'est le lien clé entre la biologie et le quantique
    rotation_angle = final_cry2_level * np.pi  # Angle max = pi (porte NOT)
    H = rotation_angle * sigmax()  # Hamiltonien pour la rotation
    
    psi0 = basis(2, 0)  # État initial |0>
    times = np.linspace(0, 1, 101) # Temps de l'opération de porte
    result = mesolve(H, psi0, times, [], [])
    final_state = result.states[-1]
    prob_1 = np.abs((basis(2, 1).dag() * final_state).full())**2

    # 3. Simulation métabolique
    # La consommation dépend du nombre d'opérations (ici, 1 porte)
    n_ops = 1
    t_metab, atp, nadph = simulate_metabolism(n_ops, sim_duration)

    return {
        "em_frequency": em_frequency,
        "final_cry2_level": final_cry2_level,
        "final_qubit_state": final_state,
        "readout_probability": prob_1,
        "metabolism": {
            "time": t_metab,
            "atp": atp,
            "nadph": nadph
        }
    }

if __name__ == '__main__':
    # Scénario de simulation
    em_freq = 80  # Hz
    results = run_multiphysics_simulation(em_freq)

    # Affichage des résultats
    print(f"--- Résultats de la simulation multiphysique (Fréquence EM = {em_freq} Hz) ---")
    print(f"Niveau d'expression final de CRY2 : {results['final_cry2_level']:.3f}")
    print(f"Probabilité de lecture de l'état |1> : {results['readout_probability']:.3f}")

    # Visualisation métabolique
    metab_results = results['metabolism']
    plt.figure(figsize=(12, 6))
    plt.plot(metab_results['time'], metab_results['atp'], label='ATP')
    plt.plot(metab_results['time'], metab_results['nadph'], label='NADPH')
    plt.title("Consommation métabolique résultante")
    plt.xlabel("Temps (s)")
    plt.ylabel("Niveau (unité arbitraire)")
    plt.legend()
    plt.grid(True)
    plt.show()
