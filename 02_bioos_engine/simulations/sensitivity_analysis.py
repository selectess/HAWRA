
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


# =============================================================================
# ANALYSE DE SENSIBILITÉ DU PROTOCOLE DE RÉGÉNÉRATION
# =============================================================================
#
# OBJECTIF :
# Identifier les étapes les plus critiques du protocole en quantifiant l'impact
# de la variation de chaque probabilité de succès sur le rendement final.
#
# MÉTHODOLOGIE :
# Analyse de sensibilité "un facteur à la fois" (OFAT). Chaque paramètre (p_i)
# est varié sur une plage définie tandis que les autres sont maintenus à leur
# valeur de base. Le rendement moyen est calculé pour chaque point de la plage.
#
# =============================================================================

def run_regeneration_simulation(n_simulations, n_explants_initial, probabilities):
    """
    Exécute une simulation de Monte Carlo du protocole de régénération.

    Args:
        n_simulations (int): Nombre de fois où l'expérience complète est simulée.
        n_explants_initial (int): Nombre d'explants de Ficus au début de chaque simulation.
        probabilities (dict): Dictionnaire des probabilités de succès pour chaque étape.

    Returns:
        list: Une liste contenant le nombre final de plantes viables pour chaque simulation.
    """
    p_transfo = probabilities['p_transfo']
    p_selection = probabilities['p_selection']
    p_callogenese = probabilities['p_callogenese']
    p_organogenese = probabilities['p_organogenese']
    p_enracinement = probabilities['p_enracinement']
    p_acclimatation = probabilities['p_acclimatation']

    final_plant_counts = []

    for _ in range(n_simulations):
        n_surviving_selection = np.random.binomial(n_explants_initial, p_transfo * p_selection)
        n_calli = np.random.binomial(n_surviving_selection, p_callogenese)
        n_shoots = np.random.binomial(n_calli, p_organogenese)
        n_rooted_plantlets = np.random.binomial(n_shoots, p_enracinement)
        n_final_plants = np.random.binomial(n_rooted_plantlets, p_acclimatation)
        final_plant_counts.append(n_final_plants)

    return final_plant_counts

def run_sensitivity_analysis(base_probabilities, n_simulations, n_explants_initial):
    """
    Exécute l'analyse de sensibilité sur les paramètres du protocole.
    """
    sensitivity_results = {}
    param_range = np.linspace(0.1, 1.0, 10) # Plage de variation de 10% à 100%

    for param_name in base_probabilities.keys():
        print(f"--- Analyse du paramètre : {param_name} ---")
        yields = []
        
        current_probabilities = base_probabilities.copy()
        for p_value in param_range:
            current_probabilities[param_name] = p_value
            
            simulation_results = run_regeneration_simulation(
                n_simulations, n_explants_initial, current_probabilities
            )
            
            mean_yield = np.mean(simulation_results) / n_explants_initial * 100
            yields.append(mean_yield)
            print(f"  {param_name} = {p_value:.2f} -> Rendement moyen = {mean_yield:.2f}%")

        sensitivity_results[param_name] = yields
    
    return sensitivity_results, param_range


def plot_sensitivity_results(sensitivity_results, param_range):
    """
    Affiche les résultats de l'analyse de sensibilité.
    """
    plt.figure(figsize=(14, 8))
    
    for param_name, yields in sensitivity_results.items():
        plt.plot(param_range, yields, marker='o', linestyle='-', label=param_name)

    plt.title("Analyse de Sensibilité des Paramètres du Protocole de Régénération")
    plt.xlabel("Probabilité de Succès de l'Étape")
    plt.ylabel("Rendement Global Moyen (%)")
    plt.legend(title="Paramètres")
    plt.grid(True)
    plt.ylim(0, max(max(y) for y in sensitivity_results.values()) * 1.1) # Ajuste l'axe Y

    output_path = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/sensitivity_analysis.png"
    plt.savefig(output_path)
    print(f"\nGraphique de l'analyse de sensibilité sauvegardé dans : {output_path}")


if __name__ == "__main__":
    # --- Paramètres de base pour la simulation ---
    N_SIMULATIONS = 2000  # Réduit pour une exécution plus rapide de l'analyse
    N_EXPLANTS_INITIAL = 500

    # --- Probabilités de base (point de fonctionnement) ---
    base_probabilities = {
        'p_transfo': 0.10,
        'p_selection': 0.60,
        'p_callogenese': 0.50,
        'p_organogenese': 0.40,
        'p_enracinement': 0.70,
        'p_acclimatation': 0.50
    }

    # --- Exécution de l'analyse de sensibilité ---
    sensitivity_data, param_range = run_sensitivity_analysis(
        base_probabilities, N_SIMULATIONS, N_EXPLANTS_INITIAL
    )

    # --- Visualisation des résultats ---
    plot_sensitivity_results(sensitivity_data, param_range)
