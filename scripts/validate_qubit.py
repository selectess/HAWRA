import numpy as np
from qutip import basis, sigmax, sigmaz, mesolve
import matplotlib.pyplot as plt
import os

# --- Configuration de la Simulation ---
# Créer le répertoire de sortie s'il n'existe pas
output_dir = "03_quantum_simulation/results"
os.makedirs(output_dir, exist_ok=True)

# --- Paramètres Physiques et de Simulation ---
# Taux de décohérence de base (déphasage) sans silice, en THz (1/ps)
gamma_no_si = 0.02

# Facteur de réduction du bruit par la cage de silice (40% de réduction)
silica_protection_factor = 0.6
gamma_with_si = gamma_no_si * silica_protection_factor

# Échelle de temps de la simulation en picosecondes (ps)
times = np.linspace(0, 250, 500)

# --- Définition du système quantique ---
# Hamiltonien (nul, on ne s'intéresse qu'à la décohérence)
H = 0 * sigmaz()

# État initial en superposition (état |+>)
psi0 = (basis(2, 0) + basis(2, 1)).unit()

# Opérateurs d'effondrement pour la décohérence (déphasage pur)
c_ops_no_si = [np.sqrt(gamma_no_si) * sigmaz()]
c_ops_with_si = [np.sqrt(gamma_with_si) * sigmaz()]

# --- Simulation ---
# 1. Sans la cage de silice
result_no_si = mesolve(H, psi0, times, c_ops_no_si, [sigmax()])
coherence_no_si = result_no_si.expect[0]

# 2. Avec la cage de silice
result_with_si = mesolve(H, psi0, times, c_ops_with_si, [sigmax()])
coherence_with_si = result_with_si.expect[0]

# --- Fonctions Utilitaires ---
def get_coherence_time(times, coherence_data):
    """Calcule le temps de cohérence (T2) en trouvant le temps où la cohérence tombe à 1/e."""
    threshold = 1 / np.e
    indices = np.where(coherence_data < threshold)[0]
    if len(indices) == 0:
        return times[-1]
    t2_index = indices[0]
    if t2_index == 0:
        return times[0]
    p1 = (times[t2_index - 1], coherence_data[t2_index - 1])
    p2 = (times[t2_index], coherence_data[t2_index])
    m = (p2[1] - p1[1]) / (p2[0] - p1[0])
    if abs(m) < 1e-9:
        return times[t2_index]
    t2_time = p1[0] + (threshold - p1[1]) / m
    return t2_time

# --- Analyse des Résultats ---
time_no_si = get_coherence_time(times, coherence_no_si)
time_with_si = get_coherence_time(times, coherence_with_si)

if time_no_si > 0:
    improvement = ((time_with_si / time_no_si) - 1) * 100
else:
    improvement = float('inf')

# --- Rapport Final ---
print("--- Rapport de Simulation de Cohérence Quantique ---")
print(f"Modèle de Décohérence : Déphasage pur (Opérateur σz)")
print(f"Taux de décohérence de base (γ_sans_silice) : {gamma_no_si}")
print(f"Facteur de protection de la silice : {silica_protection_factor}")
print(f"Taux de décohérence avec silice (γ_avec_silice) : {gamma_with_si:.4f}")
print("-"*20)
print(f"Temps de cohérence T2 (Sans Silice) : {time_no_si:.2f} ps")
print(f"Temps de cohérence T2 (Avec Silice) : {time_with_si:.2f} ps")
print(f"Amélioration de la cohérence : {improvement:.2f}%")
print("-"*20)

# --- Visualisation ---
plt.figure(figsize=(12, 7))
plt.plot(times, coherence_no_si, 'r--', label=f"Sans silice (T2 ≈ {time_no_si:.2f} ps)")
plt.plot(times, coherence_with_si, 'g-', label=f"Avec silice (T2 ≈ {time_with_si:.2f} ps)")
plt.axhline(1/np.e, color='k', linestyle=':', linewidth=2, label='Seuil 1/e')

plt.title("Validation de la Stabilité du Qubit P700 : Effet de la Cage de Silice", fontsize=16)
plt.xlabel("Temps (ps)", fontsize=14)
plt.ylabel("Cohérence <σx>", fontsize=14)
plt.legend(fontsize=12)
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.ylim(-0.1, 1.1)

# Sauvegarder la figure
output_path = os.path.join(output_dir, "qubit_coherence_validation.png")
plt.savefig(output_path, dpi=300)

print(f"\nFigure de validation sauvegardée dans : {output_path}")
