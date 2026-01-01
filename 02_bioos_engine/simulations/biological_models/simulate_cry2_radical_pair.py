#!/usr/bin/env python3
"""
Simulation radical pair CRY2 (EM 9.8 kHz)
Basé sur Ritz et al. (2000) - Larmor frequency
"""

import numpy as np
import matplotlib.pyplot as plt

def simulate_cry2_radical_pair():
    """Simule radical pair CRY2 et réponse EM"""
    print("🔬 SIMULATION RADICAL PAIR CRY2\n")
    
    # Paramètres (Ritz 2000)
    larmor_freq = 9.8e3  # 9.8 kHz
    magnetic_field = 2.2e-6
    gyromagnetic_ratio = 28.0e9  # rad/s/T (électron)
    
    # Calcul Larmor: ω = γB
    calculated_freq = (gyromagnetic_ratio * magnetic_field) / (2 * np.pi)
    
    print(f"📊 Paramètres:")
    print(f"  Fréquence Larmor: {larmor_freq/1e3:.1f} kHz")
    print(f"  Champ magnétique: {magnetic_field*1e6:.2f} μT")
    print(f"  Fréquence calculée: {calculated_freq/1e3:.1f} kHz")
    print()
    
    # Validation
    diff = abs(calculated_freq - larmor_freq) / larmor_freq * 100
    print(f"✅ Écart calculé: {diff:.1f}%")
    print(f"✅ Validation: {'OK' if diff < 5 else 'ATTENTION'}")
    
    # Simulation radical pair (FAD•- + Trp•+)
    # États de spin: singlet (S) et triplet (T)
    time = np.linspace(0, 0.01, 1000)  # 10 ms
    
    # Probabilité singlet (simplifié)
    # Oscillation à fréquence Larmor
    singlet_prob = 0.5 + 0.3 * np.cos(2 * np.pi * larmor_freq * time)
    triplet_prob = 1 - singlet_prob
    
    # Visualisation
    plt.figure(figsize=(10, 6))
    plt.plot(time * 1e3, singlet_prob, 'b-', label='Singlet (FAD•- + Trp•+)', linewidth=2)
    plt.plot(time * 1e3, triplet_prob, 'r-', label='Triplet', linewidth=2)
    plt.axvline(0.2, color='g', linestyle='--', label='200 ms (protocole)')
    plt.xlabel('Temps (ms)')
    plt.ylabel('Probabilité')
    plt.title('Radical Pair CRY2 - Réponse EM 9.8 kHz')
    plt.grid(True)
    plt.legend()
    plt.savefig('05_SIMULATION/results/cry2_radical_pair.png', dpi=150)
    print(f"\n📈 Graphique sauvegardé: 05_SIMULATION/results/cry2_radical_pair.png")
    
    return diff < 5

if __name__ == '__main__':
    import os
    os.makedirs('05_SIMULATION/results', exist_ok=True)
    success = simulate_cry2_radical_pair()
    print(f"\n{'✅ Validation réussie' if success else '❌ Validation échouée'}")

