
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# SIMULATION DU PROTOCOLE DE RÉGÉNÉRATION DE FICUS ELASTICA (HAWRA)
# =============================================================================
#
# OBJECTIF :
# Estimer le rendement global du protocole de régénération par une approche
# de Monte Carlo. La simulation modélise chaque étape critique comme un
# événement probabiliste pour prédire le nombre de plantes viables
# (HAWRA-Ficus-G0) obtenues à partir d'un nombre initial d'explants.
#
# MÉTHODOLOGIE :
# Le protocole est une chaîne de Markov où la sortie d'une étape devient
# l'entrée de la suivante. Chaque étape a une probabilité de succès
# associée, basée sur la littérature scientifique et l'expertise.
#
# PARAMÈTRES CLÉS (Probabilités de succès) :
# 1. p_transfo: Probabilité qu'une cellule d'explant intègre l'ADN-T.
# 2. p_selection: Probabilité qu'une cellule transformée survive à la sélection.
# 3. p_callogenese: Probabilité qu'un cal se forme à partir de cellules sélectionnées.
# 4. p_organogenese: Probabilité qu'un cal génère des bourgeons viables.
# 5. p_enracinement: Probabilité qu'une pousse développe un système racinaire.
# 6. p_acclimatation: Probabilité qu'une plantule survive au transfert en serre.
#
# =============================================================================

def run_regeneration_simulation(n_simulations, n_explants_initial):
    """
    Exécute une simulation de Monte Carlo du protocole de régénération.

    Args:
        n_simulations (int): Nombre de fois où l'expérience complète est simulée.
        n_explants_initial (int): Nombre d'explants de Ficus au début de chaque simulation.

    Returns:
        list: Une liste contenant le nombre final de plantes viables pour chaque simulation.
    """

    # --- Paramètres de probabilité (basés sur des estimations de la littérature) ---
    # Ces valeurs devront être affinées par des expériences pilotes.
    p_transfo = 0.10          # 10% des cellules sont transformées avec succès
    p_selection = 0.60        # 60% des cellules transformées survivent à la kanamycine
    p_callogenese = 0.50      # 50% des groupes de cellules survivants forment un cal
    p_organogenese = 0.40     # 40% des cals produisent des bourgeons
    p_enracinement = 0.70     # 70% des bourgeons développent des racines
    p_acclimatation = 0.50    # 50% des plantules survivent à l'acclimatation

    final_plant_counts = []

    for _ in range(n_simulations):
        # --- Simulation d'une seule expérience complète ---

        # Étape 1 & 2: Transformation et Sélection
        # On suppose qu'un explant peut donner naissance à plusieurs lignées cellulaires.
        # Simplification : on considère le nombre d'explants comme le nombre d'unités de départ.
        n_surviving_selection = np.random.binomial(n_explants_initial, p_transfo * p_selection)

        # Étape 3: Callogenèse
        n_calli = np.random.binomial(n_surviving_selection, p_callogenese)

        # Étape 4: Organogenèse (chaque cal peut produire plusieurs bourgeons, ici on simplifie 1 cal -> 1 pousse potentielle)
        n_shoots = np.random.binomial(n_calli, p_organogenese)

        # Étape 5: Enracinement
        n_rooted_plantlets = np.random.binomial(n_shoots, p_enracinement)

        # Étape 6: Acclimatation
        n_final_plants = np.random.binomial(n_rooted_plantlets, p_acclimatation)

        final_plant_counts.append(n_final_plants)

    return final_plant_counts

def plot_results(simulation_results, n_explants_initial):
    """
    Affiche les résultats de la simulation sous forme d'histogramme.
    """
    mean_plants = np.mean(simulation_results)
    std_plants = np.std(simulation_results)
    rendement_moyen = (mean_plants / n_explants_initial) * 100

    plt.figure(figsize=(12, 7))
    plt.hist(simulation_results, bins=range(0, max(simulation_results) + 2), alpha=0.75, edgecolor='black')
    plt.title(f"Distribution du Nombre de Plantes Viables (N_explants = {n_explants_initial})")
    plt.xlabel("Nombre de Plantes HAWRA-Ficus-G0 Viables")
    plt.ylabel("Fréquence (sur {} simulations)".format(len(simulation_results)))
    plt.axvline(mean_plants, color='r', linestyle='dashed', linewidth=2, label=f"Moyenne: {mean_plants:.2f}")
    plt.legend()

    stats_text = (
        f"--- Statistiques de la Simulation ---\n"
        f"Nombre d'explants initial: {n_explants_initial}\n"
        f"Nombre de simulations: {len(simulation_results)}\n"
        f"Nombre moyen de plantes finales: {mean_plants:.2f} ± {std_plants:.2f}\n"
        f"Rendement global moyen: {rendement_moyen:.2f}%\n\n"
        f"Conclusion : Pour {n_explants_initial} explants de départ, on peut s'attendre à obtenir \n"
        f"en moyenne {mean_plants:.2f} plantes HAWRA viables."
    )
    print(stats_text)
    
    # Sauvegarder le graphique
    output_path = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/regeneration_simulation_yield.png"
    plt.savefig(output_path)
    print(f"\nGraphique sauvegardé dans : {output_path}")


if __name__ == "__main__":
    # --- Paramètres de la simulation ---
    N_SIMULATIONS = 10000  # Nombre d'essais pour la robustesse statistique
    N_EXPLANTS_INITIAL = 500 # Nombre d'explants de départ, un nombre réaliste pour une expérience en labo

    # --- Exécution et visualisation ---
    results = run_regeneration_simulation(N_SIMULATIONS, N_EXPLANTS_INITIAL)
    plot_results(results, N_EXPLANTS_INITIAL)

