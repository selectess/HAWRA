import argparse
import argparse
import json
import os
import numpy as np
from src.hawra_simulator.simulator import Simulator

class NumpyJSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        elif isinstance(o, np.floating):
            return float(o)
        elif isinstance(o, np.ndarray):
            return o.tolist()
        elif isinstance(o, complex):
            return {'real': o.real, 'imag': o.imag}
        return super(NumpyJSONEncoder, self).default(o)

def main():
    """
    Point d'entrée principal pour le simulateur unifié HAWRA.
    """
    # Déterminer le répertoire du script pour construire des chemins absolus
    script_dir = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(description="Simulateur Unifié HAWRA-PQPE")
    parser.add_argument(
        '--config',
        type=str,
        help='Chemin vers le fichier de script .bsim.json à exécuter.'
    )
    parser.add_argument(
        '--duration',
        type=float,
        default=1e-6,
        help='Durée totale de la simulation en secondes.'
    )
    parser.add_argument(
        '--dt',
        type=float,
        default=1e-8,
        help='Pas de temps de la simulation en secondes.'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='output.json',
        help='Chemin vers le fichier de sortie pour sauvegarder le log.'
    )
    args = parser.parse_args()

    if not args.config:
        parser.error("L'argument --config est requis.")

    # Initialisation des variables
    simulator = None

    # Charger la configuration ou le script
    try:
        print(f"Chargement de la configuration depuis '{args.config}'...")
        simulator = Simulator(bsim_script=args.config)
    except FileNotFoundError:
        print(f"Erreur: Le fichier de configuration '{args.config}' n'a pas été trouvé.")
        return
    except json.JSONDecodeError:
        print(f"Erreur: Le fichier de configuration '{args.config}' n'est pas un JSON valide.")
        return

    if simulator:
        print("Initialisation du simulateur...")
        results = simulator.run()

        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=4, cls=NumpyJSONEncoder)
            print(f"Résultats de la simulation sauvegardés dans '{args.output}'.")
        except IOError as e:
            print(f"Erreur lors de la sauvegarde des résultats : {e}")

if __name__ == "__main__":
    main()
