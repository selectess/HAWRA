#!/usr/bin/env python3
"""
Simulation Épigénétique HAWRA
Simulation de l'expression génique vs température pour stabilité qubits
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

# Génome Hawra (exemple tronqué)
seq = Seq("ACGTATGGAACGAAGACGACGACGC")  # psaA tronqué
record = SeqRecord(seq, id="psaA_qubit", description="P700 chaud")

# Ajout promoteur (CaMV 35S)
prom = Seq("GAGCTCCTCCAAGAA")
full_seq = prom + seq
record.seq = full_seq

# Sauvegarde FASTA
with open("hawra.fasta", "w") as f:
    SeqIO.write(record, f, "fasta")

print("✅ Génome HAWRA sauvegardé: hawra.fasta")

# Simulation épigénétique : expression vs température
temps = np.linspace(25, 45, 21)  # 25 à 45°C
stability = 1 - 0.02 * (temps - 35) ** 2  # HSP70 garde 100% à 35°C
stability = np.clip(stability, 0, 1)  # Limite [0,1]

print("\n📊 Expression P700 selon T (°C):")
print("=" * 50)
for t, s in zip(temps, stability):
    qubits = int(1000 * s)  # 1000 qubits max
    print(f"T = {t:5.1f}°C -> {s:.2%} stabilité -> {qubits:4d} qubits actifs")

# Graphique
plt.figure(figsize=(10, 6))
plt.plot(temps, stability * 100, 'b-', linewidth=2, label='Stabilité qubits (%)')
plt.axvline(x=35, color='r', linestyle='--', label='Optimum HSP70 (35°C)')
plt.xlabel('Température (°C)', fontsize=12)
plt.ylabel('Stabilité qubits (%)', fontsize=12)
plt.title('HAWRA: Stabilité qubits vs Température (HSP70)', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.legend()
plt.ylim(0, 105)
plt.savefig('hawra_stability.png', dpi=300, bbox_inches='tight')
print("\n✅ Graphique sauvegardé: hawra_stability.png")

# Simulation cohérence quantique
coherence_time = 650 * np.exp(-(temps - 25) / 20)  # Décroissance exponentielle
coherence_time = np.clip(coherence_time, 100, 650)  # Limite [100, 650] fs

plt.figure(figsize=(10, 6))
plt.plot(temps, coherence_time, 'g-', linewidth=2, label='Temps cohérence (fs)')
plt.axhline(y=100, color='r', linestyle='--', label='Minimum fonctionnel (100 fs)')
plt.xlabel('Température (°C)', fontsize=12)
plt.ylabel('Temps cohérence (fs)', fontsize=12)
plt.title('HAWRA: Cohérence quantique P700 vs Température', fontsize=14, fontweight='bold')
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig('hawra_coherence.png', dpi=300, bbox_inches='tight')
print("✅ Graphique sauvegardé: hawra_coherence.png")

print("\n🎯 Simulation terminée!")

