
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, sigmax, sigmaz, mesolve, Bloch, Options
from matplotlib.animation import FuncAnimation

# --- Paramètres de la simulation ---
gamma = 0.1  # Taux de décohérence (déphasage)
times = np.linspace(0, 50, 200)  # Échelle de temps

# --- Définition du système quantique ---
# Hamiltonien (nul, on ne s'intéresse qu'à la décohérence)
H = 0 * sigmaz()

# État initial : superposition |+> = (|0> + |1>)/sqrt(2)
psi0 = (basis(2, 0) + basis(2, 1)).unit()

# Opérateurs d'effondrement pour la décohérence (déphasage pur)
# c_ops simule l'interaction avec l'environnement
c_ops = [np.sqrt(gamma) * sigmaz()]

# --- Exécution de la simulation ---
# On veut observer la décroissance de la cohérence et stocker les états pour l'animation.
options = Options(store_states=True)
result = mesolve(H, psi0, times, c_ops, [sigmax()], options=options)

# --- Création de l'animation de la sphère de Bloch ---
fig = plt.figure(figsize=(8, 8))
b = Bloch(fig=fig)

def animate(i):
    b.clear()
    b.add_states(result.states[i])
    b.make_sphere()
    return b.fig.axes[0].artists

animation = FuncAnimation(fig, animate, frames=len(times), blit=True)

# Sauvegarde de l'animation en GIF
animation.save('bloch_sphere_decoherence.gif', writer='imagemagick', fps=20, dpi=100)
print("Animation de la sphère de Bloch sauvegardée dans 'bloch_sphere_decoherence.gif'")

# --- Analyse et Visualisation (graphique 2D) ---
# La valeur attendue de sigmax() est une mesure de la cohérence de phase
coherence = result.expect[0]

# Calcul théorique pour comparaison
theoretical_coherence = np.exp(-2 * gamma * times)

def get_coherence_time(times, coherence_data):
    """Calcule le temps de cohérence (T2) où la cohérence tombe à 1/e."""
    threshold = 1 / np.e
    # Trouve l'index où la cohérence passe sous le seuil
    indices = np.where(coherence_data < threshold)[0]
    if len(indices) > 0:
        # Prend le premier point de passage
        t2_index = indices[0]
        # Interpolation linéaire simple pour une meilleure précision
        if t2_index > 0:
            p1 = (times[t2_index - 1], coherence_data[t2_index - 1])
            p2 = (times[t2_index], coherence_data[t2_index])
            # (y - y1) = m * (x - x1) => x = x1 + (y - y1) / m
            m = (p2[1] - p1[1]) / (p2[0] - p1[0])
            t2_time = p1[0] + (threshold - p1[1]) / m
            return t2_time
        else:
            return times[t2_index]
    return 0 # Retourne 0 si la cohérence ne tombe jamais sous le seuil

# Calcul du temps de cohérence T2
t2_calculated = get_coherence_time(times, coherence)
# Pour un déphasage pur avec sigmaz, la décroissance de <sigmax> est en exp(-2*gamma*t)
t2_theoretical = 1 / (2 * gamma)

print("Simulation terminée.")
print(f"Valeur de cohérence initiale (attendue ~1.0): {coherence[0]:.4f}")
print(f"Valeur de cohérence finale (attendue ~0.0): {coherence[-1]:.4f}")
print("-" * 30)
print(f"Temps de cohérence T2 calculé : {t2_calculated:.4f}")
print(f"Temps de cohérence T2 théorique (1 / (2*gamma)): {t2_theoretical:.4f}")
print("-" * 30)


# --- Visualisation 2D ---
plt.figure(figsize=(10, 6))
plt.plot(times, coherence, label='Coherence simulée (<sigmax>)')
plt.plot(times, theoretical_coherence, 'r--', label='Décroissance théorique (exp(-2*gamma*t))')
plt.axhline(1/np.e, color='g', linestyle=':', label='Seuil 1/e')
plt.axvline(t2_calculated, color='purple', linestyle='--', label=f'T2 calculé = {t2_calculated:.2f}')
plt.title("Débogage de la Décohérence avec QuTiP (Déphasage pur)")
plt.xlabel("Temps")
plt.ylabel("Cohérence (<sigmax>)")
plt.legend()
plt.grid(True)
plt.savefig("debug_coherence.png")
print("Graphique de débogage sauvegardé dans 'debug_coherence.png'")
plt.show()
