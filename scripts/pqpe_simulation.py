# %% [markdown]
# # Simulation Numérique du Protocole PQPE HAWRA
#
# **Projet :** HAWRA
# **Standard :** Zéro Tolérance - Haute Exigence Scientifique
#
# ## 1. Objectif
#
# Ce script simule les phases opérationnelles (II, III, IV) du protocole PQPE en utilisant des unités physiques réelles et des modèles validés, conformément au Protocole Zéro Tolérance.
#
# - **Phase II :** Calibration des portes quantiques.
# - **Phase III :** Exécution d'un algorithme simple.
# - **Phase IV :** Simulation de la dérive biologique et de la re-calibration.

# %%
# --- Import des bibliothèques ---
import numpy as np
import matplotlib.pyplot as plt
from qutip import basis, sigmax, sigmay, sigmaz, sigmam, mesolve, fidelity, Qobj, tensor
import os

print("Bibliothèques importées.")

def define_physical_system():
    """Définit les paramètres physiques du bio-qubit P700 et les états clés."""
    # Crée le répertoire de résultats s'il n'existe pas
    results_dir = '05_data/results/pqpe_simulation'
    os.makedirs(results_dir, exist_ok=True)

    # Temps de cohérence T2 (en nanosecondes)
    T2_time = 1000.0  # ns

    # Taux de déphasage correspondant (opérateur sigma-z pour le déphasage)
    # gamma_phi = 1 / T2. L'opérateur est sqrt(gamma_phi) * sigmaz
    gamma_phi = 1.0 / T2_time
    
    # Taux de relaxation T1 (supposé plus long que T2)
    T1_time = 2000.0 # ns
    gamma_relax = 1.0 / T1_time

    # Opérateurs d'effondrement : relaxation (T1) et déphasage (T2)
    c_ops = [
        np.sqrt(gamma_relax) * sigmam(),  # Relaxation T1
        np.sqrt(gamma_phi) * sigmaz()      # Déphasage T2
    ]

    # États de base
    initial_state = basis(2, 0)  # État |0>
    target_state = basis(2, 1)   # État |1> (cible pour une porte NOT)
    
    # Retourner les paramètres essentiels
    # Le gamma retourné est le déphasage, car c'est souvent le plus critique
    return initial_state, target_state, c_ops, gamma_phi


# %% [markdown]
# ## 2. Définition du Système Physique (Unités Réelles)
#
# Nous définissons ici le modèle physique du bio-qubit P700.

# %%

# %% [markdown]
# ## 3. Phase II : Calibration des Portes Quantiques
#
# ### 3.1. Définition des impulsions de contrôle
#
# Nous simulons des portes quantiques en appliquant des impulsions micro-ondes résonnantes. Une impulsion avec une enveloppe gaussienne est un choix réaliste.

# %%

def calibrate_gate(initial_state, target_state, c_ops):
    """
    Phase II: Calibre la porte quantique (NOT gate) en trouvant les paramètres d'impulsion optimaux.
    """
    print("--- Début de la Calibration (Phase II) ---")
    
    # Espace de recherche pour les paramètres d'impulsion (durée et amplitude)
    # Corrigé pour inclure la valeur théorique (A*t = 0.5 pour H=pi*A*sx)
    pulse_durations = np.linspace(1, 20, 20)  # ns
    pulse_amplitudes = np.linspace(0.01, 0.15, 20) # GHz

    best_fidelity = 0
    best_params = {}

    for duration in pulse_durations:
        for amplitude in pulse_amplitudes:
            # Le Hamiltonien pour une impulsion de Rabi. H = (Omega/2) * sigma_x
            # Omega (Rabi freq) = 2*pi * A. Donc H = pi * A * sigma_x
            # L'argument de l'exponentielle dans mesolve est -i*H*t.
            # Pour une porte NOT (rotation de pi), on veut que l'angle soit pi.
            # L'opérateur de rotation est exp(-i * (theta/2) * sigma_x).
            # Ici, l'angle est 2 * (pi * A * t). On veut que ce soit pi.
            # Donc 2 * pi * A * t = pi => A * t = 0.5
            H = np.pi * amplitude * sigmax() # Hamiltonien en unités de rad/ns
            times = np.linspace(0, duration, 100)
            
            result = mesolve(H, initial_state, times, c_ops, [])
            
            current_fidelity = fidelity(result.states[-1], target_state)
            
            if current_fidelity > best_fidelity:
                best_fidelity = current_fidelity
                best_params = {'duration': duration, 'amplitude': amplitude}

    print("Calibration terminée.")
    print(f"  - Meilleure fidélité atteinte : {best_fidelity:.4f}")
    print(f"  - Durée d'impulsion optimale : {best_params.get('duration', 0):.2f} ns")
    print(f"  - Amplitude d'impulsion optimale : {best_params.get('amplitude', 0):.4f} GHz")
    
    if best_fidelity < 0.95:
        print("ALERTE ZT: Fidélité de porte insuffisante après calibration (< 95%).")
    else:
        print(f"Calibration réussie avec une fidélité de {best_fidelity:.4f}.")

    return best_params

def apply_gate(params, initial_state, c_ops):
    """Applique la porte avec les paramètres donnés."""
    H = np.pi * params['amplitude'] * sigmax() # Hamiltonien en unités de rad/ns
    times = np.linspace(0, params['duration'], 100)
    result = mesolve(H, initial_state, times, c_ops, [])
    return result.states[-1]

def get_gate_fidelity(params, initial_state, target_state, c_ops):
    """Calcule la fidélité de la porte pour des paramètres donnés."""
    if not params:
        return 0
    final_state = apply_gate(params, initial_state, c_ops)
    return fidelity(final_state, target_state)

# --- Exécution principale de la simulation ---
if __name__ == "__main__":
    # Définition du système physique (Phase I implicite)
    initial_state, target_state, c_ops, gamma = define_physical_system()
    print("Système physique défini.")

    # Phase II: Calibration
    optimal_params = calibrate_gate(initial_state, target_state, c_ops)

    # Phase III: Opération - Validation de la porte calibrée
    print("\n--- Validation de la porte calibrée (Phase III) ---")
    op_fidelity = get_gate_fidelity(optimal_params, initial_state, target_state, c_ops)
    print(f"Exécution de la porte X calibrée sur l'état |0>.")
    print(f"  - Fidélité de l'état final par rapport à |1> : {op_fidelity:.4f}")
    if op_fidelity < 0.95:
        print("ALERTE ZT: La porte calibrée ne fonctionne pas comme prévu.")

    # Phase IV: Maintenance - Simulation de la dérive et re-calibration
    print("\n--- Simulation de la dérive et re-calibration (Phase IV) ---")
    
    # Simuler une dérive (e.g., augmentation du bruit)
    gamma_drifted = gamma * 2.5 
    c_ops_drifted = [np.sqrt(gamma_drifted) * sigmam()]
    print(f"Dérive biologique simulée : le taux de décohérence a augmenté à {gamma_drifted:.4f}.")

    # Vérifier la performance avec les anciens paramètres optimaux
    fidelity_drifted = get_gate_fidelity(optimal_params, initial_state, target_state, c_ops_drifted)
    print("Performance avec dérive (anciens paramètres):")
    print(f"  - Nouvelle fidélité de la porte X : {fidelity_drifted:.4f} (contre {op_fidelity:.4f} initialement)")

    # Lancer la re-calibration sur le système dérivé
    print("Lancement de la re-calibration...")
    recalibrated_params = calibrate_gate(initial_state, target_state, c_ops_drifted)

    # Vérifier la performance après re-calibration
    fidelity_recalibrated = get_gate_fidelity(recalibrated_params, initial_state, target_state, c_ops_drifted)
    print("\nPerformance après re-calibration:")
    print(f"  - Fidélité restaurée : {fidelity_recalibrated:.4f}")
    if recalibrated_params:
        print(f"  - Nouveaux paramètres optimaux : {recalibrated_params['duration']:.2f} ns, {recalibrated_params['amplitude']:.4f} GHz")

    print("\n--- Fin de la simulation du protocole PQPE ---")