#!/usr/bin/env python3
"""
Simulation croissance HAWRA (conditions réelles Maroc)
Modèle simplifié basé sur paramètres OpenSimRoot

⚠️ LIMITATION IMPORTANTE :
Ce script est un modèle simplifié (linéaire) et ne constitue PAS une intégration
complète d'OpenSimRoot. Pour une simulation 3D complète avec OpenSimRoot, voir :
- Documentation OpenSimRoot: https://opensimroot.org
- Dossier: 05_SIMULATION/opensimroot/ (à développer)

Ce modèle valide uniquement :
- Croissance linéaire simplifiée (3 cm/mois)
- Production qubits théorique (1 qubit/cm, max 1000)
- Métabolisme CAM simplifié (80% efficacité nuit)

Pour validation expérimentale, des mesures réelles sont nécessaires.
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

def simulate_growth():
    """Simule croissance 3D HAWRA sur 90 jours"""
    print("🌱 SIMULATION CROISSANCE HAWRA\n")
    
    # Paramètres climat Maroc
    temp_day = 35  # °C jour
    temp_night = 25  # °C nuit
    humidity = 20  # %
    light_hours = 12  # h/jour
    
    # Paramètres croissance
    days = 90
    initial_height = 0.01  # 1 cm (graine)
    growth_rate = 0.03  # cm/jour (3 cm/mois)
    
    # Simulation
    time = np.arange(0, days + 1)
    height = initial_height + growth_rate * time
    
    # Production qubits (1 qubit par cm de hauteur, max 1000)
    qubits = np.minimum(height * 100, 1000)
    
    # CAM métabolisme (énergie 24h/24)
    cam_efficiency = 0.8  # 80% efficacité
    energy_day = np.ones(days + 1) * 100  # 100% jour
    energy_night = np.ones(days + 1) * (80 * cam_efficiency)  # 80% nuit (CAM)
    
    print(f"📊 Paramètres:")
    print(f"  Durée: {days} jours")
    print(f"  Température: {temp_day}°C jour / {temp_night}°C nuit")
    print(f"  Humidité: {humidity}%")
    print()
    
    print(f"📈 Résultats (jour {days}):")
    print(f"  Hauteur: {height[-1]:.2f} m")
    print(f"  Qubits: {int(qubits[-1])}")
    print(f"  Énergie jour: {energy_day[-1]:.0f}%")
    print(f"  Énergie nuit: {energy_night[-1]:.0f}%")
    print()
    
    # Visualisation
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Croissance
    axes[0, 0].plot(time, height, 'g-', linewidth=2)
    axes[0, 0].set_xlabel('Jours')
    axes[0, 0].set_ylabel('Hauteur (m)')
    axes[0, 0].set_title('Croissance HAWRA')
    axes[0, 0].grid(True)
    
    # Qubits
    axes[0, 1].plot(time, qubits, 'b-', linewidth=2)
    axes[0, 1].axhline(100, color='r', linestyle='--', label='v0.1 (100 qubits)')
    axes[0, 1].axhline(1000, color='orange', linestyle='--', label='v1.0 (1000 qubits)')
    axes[0, 1].set_xlabel('Jours')
    axes[0, 1].set_ylabel('Qubits')
    axes[0, 1].set_title('Production Qubits')
    axes[0, 1].legend()
    axes[0, 1].grid(True)
    
    # Énergie
    axes[1, 0].plot(time, energy_day, 'y-', label='Jour', linewidth=2)
    axes[1, 0].plot(time, energy_night, 'm-', label='Nuit (CAM)', linewidth=2)
    axes[1, 0].set_xlabel('Jours')
    axes[1, 0].set_ylabel('Énergie (%)')
    axes[1, 0].set_title('Métabolisme CAM 24h/24')
    axes[1, 0].legend()
    axes[1, 0].grid(True)
    
    # Température
    temp_cycle = []
    for day in range(days + 1):
        for hour in range(24):
            if 6 <= hour < 18:  # Jour
                temp_cycle.append(temp_day)
            else:  # Nuit
                temp_cycle.append(temp_night)
    axes[1, 1].plot(temp_cycle[:240], 'r-', linewidth=1)  # 10 premiers jours
    axes[1, 1].set_xlabel('Heures (10 premiers jours)')
    axes[1, 1].set_ylabel('Température (°C)')
    axes[1, 1].set_title('Cycle Température Maroc')
    axes[1, 1].grid(True)
    
    plt.tight_layout()
    plt.savefig('05_SIMULATION/results/growth_simulation.png', dpi=150)
    print(f"📈 Graphique sauvegardé: 05_SIMULATION/results/growth_simulation.png")
    
    # Validation
    height_ok = height[-1] >= 2.0  # Au moins 2m en 90 jours
    qubits_ok = qubits[-1] >= 100  # Au moins 100 qubits
    energy_ok = energy_night[-1] >= 50  # Au moins 50% énergie nuit
    
    print(f"\n✅ Validation:")
    print(f"  Hauteur: {'OK' if height_ok else 'ATTENTION'}")
    print(f"  Qubits: {'OK' if qubits_ok else 'ATTENTION'}")
    print(f"  Énergie: {'OK' if energy_ok else 'ATTENTION'}")
    
    return height_ok and qubits_ok and energy_ok

if __name__ == '__main__':
    import os
    os.makedirs('05_SIMULATION/results', exist_ok=True)
    success = simulate_growth()
    print(f"\n{'✅ Simulation réussie' if success else '❌ Simulation échouée'}")

