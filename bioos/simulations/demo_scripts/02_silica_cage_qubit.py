import numpy as np

# Sans silice : 650 fs
t = np.linspace(0, 2e-12, 2000)
coherence_no_si = np.exp(-t / 650e-15)

# Avec cage silice : +20%
coherence_si = np.exp(-t / (650e-15 * 1.2))

# Sauvegarder en CSV
data = np.vstack((t, coherence_no_si, coherence_si)).T
np.savetxt("silica_cage.csv", data, delimiter=",", header="time,no_si,si")

print("Données de cohérence générées dans silica_cage.csv")
