
import json
import matplotlib.pyplot as plt
import os

def visualize_electrochemical(json_file, output_image):
    """
    Génère une visualisation des aspects électrochimiques à partir d'un fichier de simulation JSON.
    """
    with open(json_file, 'r') as f:
        data = json.load(f)

    time = [d['time'] for d in data]
    light_intensity = [d['light_intensity'] for d in data]
    p700_concentration = [d['p700_concentration'] for d in data]
    luc_green_output = [d['luc_green_output'] for d in data]

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 18), sharex=True)

    # Graphe 1: Signal de contrôle électrique (LED)
    ax1.plot(time, light_intensity, color='blue', label='Intensité lumineuse (LED)')
    ax1.set_ylabel('Intensité lumineuse (unité arbitraire)')
    ax1.set_title('Signal de Contrôle Électrique')
    ax1.legend()
    ax1.grid(True)

    # Graphe 2: État quantique (P700)
    ax2.plot(time, p700_concentration, color='red', label='Concentration de P700')
    ax2.set_ylabel('Concentration (unité arbitraire)')
    ax2.set_title('État Quantique du Système')
    ax2.legend()
    ax2.grid(True)

    # Graphe 3: Réponse électrochimique (Luciférase)
    ax3.plot(time, luc_green_output, color='green', label='Sortie de la Luciférase')
    ax3.set_xlabel('Temps (unité arbitraire)')
    ax3.set_ylabel('Luminescence (unité arbitraire)')
    ax3.set_title('Réponse Électrochimique')
    ax3.legend()
    ax3.grid(True)

    plt.suptitle('Visualisation Électrique, Quantique et Électrochimique du Système HAWRA', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Créer le répertoire de sortie s'il n'existe pas
    output_dir = os.path.dirname(output_image)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plt.savefig(output_image, dpi=300)
    print(f"Graphique électrochimique sauvegardé dans {output_image}")

if __name__ == "__main__":
    visualize_electrochemical(
        "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/multiphysics_simulation/multiphysics_simulation_v2.json",
        "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/electrochemical_visualization/electrochemical_response.png"
    )
