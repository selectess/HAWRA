
import random
import time

# --- Paramètres de la Simulation ---
# Probabilités de succès pour chaque étape critique du protocole SOP.
# Ces valeurs sont basées sur la littérature et l'expertise du domaine.
PROBABILITIES = {
    "explant_survival": 0.95,          # Survie de l'explant après prélèvement
    "transformation_efficiency": 0.40, # Efficacité de la transformation par Agrobacterium
    "selection_survival": 0.25,        # Survie à la sélection par antibiotique/herbicide
    "callus_formation": 0.70,          # Formation de cals à partir des cellules transformées
    "shoot_induction": 0.30,           # Induction de pousses à partir des cals
    "root_formation": 0.50,            # Développement de racines sur les pousses
    "acclimatization_survival": 0.60   # Survie de la plante lors du passage en terre
}

def run_procedural_simulation():
    """
    Exécute une simulation impérative du protocole de régénération pour un seul explant.
    Affiche le résultat de chaque étape.
    """
    print("===========================================================")
    print("=== Lancement de la Simulation Procédurale du SOP Ficus ===")
    print("=== Suivi du parcours d'un seul explant...              ===")
    print("===========================================================")
    time.sleep(1)

    # Étape 1: Prélèvement de l'explant
    print("\n[ÉTAPE 1/7] Prélèvement et préparation de l'explant...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["explant_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu à la stérilisation/préparation.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'explant est viable.")

    # Étape 2: Transformation génétique
    print("\n[ÉTAPE 2/7] Co-culture avec Agrobacterium pour transformation...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["transformation_efficiency"]:
        print("  -> RÉSULTAT: ÉCHEC. Le transfert du gène HAWRA a échoué.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. L'ADN-T a été transféré avec succès.")

    # Étape 3: Sélection
    print("\n[ÉTAPE 3/7] Sélection des cellules transformées sur milieu sélectif...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["selection_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. L'explant n'a pas survécu à l'agent de sélection.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. Les cellules transformées survivent et prolifèrent.")

    # Étape 4: Formation du cal
    print("\n[ÉTAPE 4/7] Induction de la formation de cals...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["callus_formation"]:
        print("  -> RÉSULTAT: ÉCHEC. Les cellules n'ont pas formé de cal embryogène.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. Un cal viable s'est développé.")

    # Étape 5: Induction des pousses
    print("\n[ÉTAPE 5/7] Induction des pousses à partir du cal...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["shoot_induction"]:
        print("  -> RÉSULTAT: ÉCHEC. Le cal n'a pas produit de pousses.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. Une ou plusieurs pousses se sont formées.")

    # Étape 6: Enracinement
    print("\n[ÉTAPE 6/7] Développement du système racinaire...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["root_formation"]:
        print("  -> RÉSULTAT: ÉCHEC. Les pousses n'ont pas développé de racines.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    print("  -> RÉSULTAT: SUCCÈS. La plantule a des racines fonctionnelles.")

    # Étape 7: Acclimatation
    print("\n[ÉTAPE 7/7] Acclimatation de la plantule en conditions ex vitro...")
    time.sleep(0.5)
    if random.random() > PROBABILITIES["acclimatization_survival"]:
        print("  -> RÉSULTAT: ÉCHEC. La plante n'a pas survécu au transfert en terre.")
        print("\n--- SIMULATION TERMINÉE ---")
        return False
    
    print("  -> RÉSULTAT: SUCCÈS. La plante est acclimatée et viable !")
    print("\n=================================================================")
    print("=== 🎉 FÉLICITATIONS ! Une plante HAWRA-Ficus-G0 a été créée ! ===")
    print("=================================================================")
    return True

if __name__ == "__main__":
    run_procedural_simulation()
