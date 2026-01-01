
import numpy as np

def simulate_cry2_expression(em_signal_frequency, duration, time_step=0.1):
    """
    Simule l'expression du gène CRY2 en réponse à un signal électromagnétique.

    Args:
        em_signal_frequency (float): Fréquence du signal EM en Hz.
        duration (float): Durée de la simulation en secondes.
        time_step (float): Pas de temps pour la simulation en secondes.

    Returns:
        tuple: Un tuple contenant les tableaux de temps et de niveau d'expression.
    """
    time = np.arange(0, duration, time_step)
    # Modèle simplifié : l'expression est proportionnelle à la fréquence du signal
    # avec une réponse de type sigmoïde pour représenter la saturation.
    # Les paramètres de la sigmoïde sont choisis pour illustrer le concept.
    k = 0.1  # Raideur de la sigmoïde
    midpoint = 50  # Fréquence à mi-expression
    expression_level = 1 / (1 + np.exp(-k * (em_signal_frequency - midpoint)))
    
    # Simuler une dynamique temporelle simple (atteinte du plateau)
    expression_over_time = expression_level * (1 - np.exp(-time / (duration / 5)))

    return time, expression_over_time

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # Exemple de simulation
    freq = 70  # Hz
    sim_duration = 10  # secondes
    t, expression = simulate_cry2_expression(freq, sim_duration)

    # Visualisation
    plt.figure(figsize=(10, 6))
    plt.plot(t, expression)
    plt.title(f"Simulation de l'expression de CRY2 (Fréquence EM = {freq} Hz)")
    plt.xlabel("Temps (s)")
    plt.ylabel("Niveau d'expression (normalisé)")
    plt.grid(True)
    plt.ylim(0, 1)
    plt.show()
