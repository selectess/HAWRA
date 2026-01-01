
import numpy as np

def simulate_metabolism(quantum_operations_count, duration, time_step=0.1, initial_atp=100, initial_nadph=100):
    """
    Simule la consommation d'ATP et de NADPH en fonction de l'activité quantique.

    Args:
        quantum_operations_count (int): Nombre d'opérations quantiques effectuées.
        duration (float): Durée de la simulation en secondes.
        time_step (float): Pas de temps pour la simulation en secondes.
        initial_atp (float): Niveau initial d'ATP.
        initial_nadph (float): Niveau initial de NADPH.

    Returns:
        tuple: Un tuple contenant les tableaux de temps, de niveaux d'ATP et de niveaux de NADPH.
    """
    time = np.arange(0, duration, time_step)
    
    # Modèle de consommation : proportionnelle au nombre d'opérations quantiques
    atp_consumption_rate = 0.1 * quantum_operations_count
    nadph_consumption_rate = 0.05 * quantum_operations_count
    
    atp_levels = initial_atp - atp_consumption_rate * time
    nadph_levels = initial_nadph - nadph_consumption_rate * time
    
    # S'assurer que les niveaux ne deviennent pas négatifs
    atp_levels[atp_levels < 0] = 0
    nadph_levels[nadph_levels < 0] = 0

    return time, atp_levels, nadph_levels

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # Exemple de simulation
    n_ops = 50  # Nombre d'opérations quantiques
    sim_duration = 20  # secondes
    t, atp, nadph = simulate_metabolism(n_ops, sim_duration)

    # Visualisation
    plt.figure(figsize=(12, 6))
    plt.plot(t, atp, label='ATP')
    plt.plot(t, nadph, label='NADPH')
    plt.title(f"Simulation de la consommation métabolique (Opérations quantiques = {n_ops})")
    plt.xlabel("Temps (s)")
    plt.ylabel("Niveau (unité arbitraire)")
    plt.legend()
    plt.grid(True)
    plt.show()
