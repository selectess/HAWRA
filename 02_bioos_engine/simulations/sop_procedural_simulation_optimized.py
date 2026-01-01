
import random
import time

# --- Paramètres de la Simulation Optimisée ---
# Probabilités de succès pour chaque étape, avec une efficacité de transformation améliorée.
PROBABILITIES = {
    "explant_survival": 0.95,          # Survie de l'explant après prélèvement
    "transformation_efficiency": 0.65, # Efficacité de transformation AMÉLIORÉE (de 0.40 à 0.65)
    "selection_survival": 0.25,        # Survie à la sélection par antibiotique/herbicide
    "callus_formation": 0.70,          # Formation de cals à partir des cellules transformées
    "shoot_induction": 0.30,           # Induction de pousses à partir des cals
    "root_formation": 0.50,            # Développement de racines sur les pousses
    "acclimatization_survival": 0.60   # Survie de la plante lors du passage en terre
}

def run_procedural_simulation_optimized():
    """
    Exécute une simulation impérative du protocole de régénération OPTIMISÉ.
    """
    print("=====================================================================")
    print("=== Lancement de la Simulation du Protocole OPTIMISÉ (200µM AS) ===")
    print("=== Suivi du parcours d'un seul explant...                       ===")
    print("=====================================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")

    # Étape 2: Transformation génétique (OPTIMISÉE)
    print("\n[ÉTAPE 2/7] Co-culture avec Agrobacterium (Protocole Optimisé, 200µM Acetosyringone)...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["transformation_efficiency"]:
        print("  -> RÉSULTAT: ÉCHEC. Le transfert du gène HAWRA a échoué.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'ADN-T a été transféré avec une meilleure efficacité.")

    # Les étapes suivantes restent identiques
    print("\n[ÉTAPE 3/7] Sélection des cellules transformées...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["selection_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu à la sélection.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. Les cellules transformées survivent.")

    print("\n[ÉTAPE 4/7] Induction de la formation de cals...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["callus_formation"]:
        print("  -> RÉSULTAT: ÉCHEC. Pas de formation de cal.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. Un cal viable s'est développé.")

    print("\n[ÉTAPE 5/7] Induction des pousses...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["shoot_induction"]:
        print("  -> RÉSULTAT: ÉCHEC. Pas de production de pousses.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. Des pousses se sont formées.")

    print("\n[ÉTAPE 6/7] Développement du système racinaire...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["root_formation"]:
        print("  -> RÉSULTAT: ÉCHEC. Pas de développement de racines.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. La plantule a des racines.")

    print("\n[ÉTAPE 7/7] Acclimatation ex vitro...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["acclimatization_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. La plante n'a pas survécu au transfert.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    
    print("  -> RÉSULTAT: SUCCÈS. La plante est acclimatée et viable !")
    print("\n=================================================================")
    print("=== 🎉 FÉLICITATIONS ! Une plante HAWRA-Ficus-G0 a été créée ! ===")
    print("=================================================================")
    return True

if __name__ == "__main__":
    run_procedural_simulation_optimized()
