
import json
import matplotlib.pyplot as plt
import numpy as np
from qutip import Bloch, basis, sigmax, sigmay, sigmaz, expect
from PIL import Image
import os

def plot_bloch_sphere(states, title, output_filename):
    """Génère une animation GIF de la sphère de Bloch pour une série d'états quantiques."""
    filenames = []
    for i, state_vector in enumerate(states):
        b = Bloch()
        # S'assurer que state_vector est une liste de nombres complexes
        state_vector_complex = [complex(s['real'], s['imag']) for s in state_vector]
        
        # Normaliser le vecteur d'état manuellement
        state_vector_np = np.array(state_vector_complex)
        norm = np.linalg.norm(state_vector_np)

        if norm < 1e-9:
            vec = [0, 0, 0]
        else:
            normalized_vector = state_vector_np / norm
            state = normalized_vector[0] * basis(2, 0) + normalized_vector[1] * basis(2, 1)
            vec = [np.real(expect(sigmax(), state)),
                   np.real(expect(sigmay(), state)),
                   np.real(expect(sigmaz(), state))]
        
        b.add_vectors(vec)
        b.make_sphere()
        plt.title(f"{title} (t={i})")
        filename = f"temp_{i}.png"
        b.save(filename)
        filenames.append(filename)
        plt.close()

    # Créer le GIF avec Pillow
    images = [Image.open(filename) for filename in filenames]
    images[0].save(output_filename,
                   save_all=True, append_images=images[1:], optimize=False, duration=100, loop=0)

    # Supprimer les fichiers temporaires
    for filename in filenames:
        os.remove(filename)

def analyze_pqpe(data):
    """Analyse le processus PQPE à partir des données de simulation."""
    if "quantum_history" in data and data["quantum_history"]:
        print("Analyse de l'historique quantique...")
        # Extraire l'historique des états du qubit 'P700_state'
        p700_states = [entry['state_vector'] for entry in data["quantum_history"]]
        
        # Générer l'animation de la sphère de Bloch
        plot_bloch_sphere(p700_states, "PQPE Process for P700 Qubit", "pqpe_bloch_sphere.gif")
        print("Animation de la sphère de Bloch générée : pqpe_bloch_sphere.gif")
    else:
        print("Aucune donnée d'historique quantique trouvée.")

if __name__ == "__main__":
    try:
        with open("output.json", "r") as f:
            simulation_data = json.load(f)
        analyze_pqpe(simulation_data)
    except FileNotFoundError:
        print("Erreur : Le fichier output.json n'a pas été trouvé. Exécutez d'abord la simulation.")
    except json.JSONDecodeError:
        print("Erreur : Impossible de décoder le fichier JSON. Le format est peut-être incorrect.")
