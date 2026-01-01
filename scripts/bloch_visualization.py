import qutip as qt
import matplotlib.pyplot as plt
from qutip import Bloch

# Simulation simple de sphère de Bloch pour état quantique (basé sur résultats PAQPE)
b = Bloch()
b.add_states(qt.basis(2, 0))  # État initial |0>
b.add_states(qt.basis(2, 0) + qt.basis(2, 1))  # Après Hadamard (superposition)
b.save('/Users/mehdiwhb/Desktop/HAWRA/bloch_sphere_paqpe.png')  # Sauvegarde récente style quantum biology

# Courbe d'intensité lumineuse vs. cohérence (dynamique)
intensities = [500, 1000, 1500]
coherences = [0.8, 0.9, 1.0]  # Exemples basés sur simulation
plt.plot(intensities, coherences)
plt.xlabel('Intensité Lumineuse')
plt.ylabel('Facteur de Cohérence')
plt.title('Robustesse à Stimuli Dynamiques - HAWRA PAQPE')
plt.savefig('/Users/mehdiwhb/Desktop/HAWRA/coherence_plot.png')
