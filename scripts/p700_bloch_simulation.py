import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, sigmax, sigmaz, mesolve, Bloch, Options
from matplotlib.animation import FuncAnimation

# --- Paramètres de la simulation pour le bio-qubit P700 ---
# Valeur hypothétique basée sur la littérature pour des systèmes similaires.
# Un temps de cohérence T2 de 1 microseconde est optimiste mais plausible.
T2_p700 = 1000  # en nanosecondes
gamma = 1 / (2 * T2_p700)  # Taux de déphasage

# Durée de la simulation
t_max = 2 * T2_p700
times = np.linspace(0, t_max, 200)

# --- Définition du système quantique ---
# Hamiltonien (pas d'évolution propre, seulement la décohérence)
H = 0 * sigmaz()

# État initial : superposition |+>
psi0 = (basis(2, 0) + basis(2, 1)).unit()

# Opérateurs de décohérence (effondrement)
c_ops = [np.sqrt(gamma) * sigmaz()]

# --- Exécution de la simulation ---
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

animation.save('p700_bloch_sphere.gif', writer='imagemagick', fps=20, dpi=100)
print("Animation de la sphère de Bloch pour P700 sauvegardée dans 'p700_bloch_sphere.gif'")

# --- Analyse et Visualisation ---
coherence = result.expect[0]
theoretical_coherence = np.exp(-2 * gamma * times)

plt.figure(figsize=(10, 6))
plt.plot(times, coherence, label='Cohérence simulée P700 (<sigmax>)')
plt.plot(times, theoretical_coherence, 'r--', label=f'Décroissance théorique (T2={T2_p700} ns)')
plt.title("Simulation de la Décohérence du Bio-Qubit P700")
plt.xlabel("Temps (ns)")
plt.ylabel("Cohérence (<σx>)")
plt.legend()
plt.grid(True)
plt.savefig('p700_coherence_decay.png')
print("Graphique de la décroissance de cohérence pour P700 sauvegardé dans 'p700_coherence_decay.png'")
