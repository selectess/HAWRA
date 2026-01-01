import os
import sys

# Correction de la logique d'importation :
# Ajout des chemins au sys.path AVANT les tentatives d'importation.
script_dir = os.path.dirname(os.path.abspath(__file__))
# Chemin vers la racine du projet (HAWRA)
project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))

# Unified simulator source path
simulator_src = os.path.join(project_root, '03_unified_simulator', 'src')

if simulator_src not in sys.path:
    sys.path.append(simulator_src)
if project_root not in sys.path:
    sys.path.append(project_root)


import time
import json

try:
    from hawra_simulator.simulator import Simulator
except ImportError:
    # Fallback pour le linter ou cas particuliers
    sys.path.append(os.path.join(project_root, '03_unified_simulator', 'src'))
    from hawra_simulator.simulator import Simulator
import subprocess
import requests

class BioOS:
    """
    Le système d'exploitation biologique (BioOS) pour la gestion d'une
    Phyto-synthetic Quantum Processing Entity (PQPE).
    """
    def __init__(self, config: dict):
        """
        Initialise le BioOS avec une configuration.

        Args:
            config (dict): Dictionnaire de configuration pour le simulateur.
        """
        print(f"[BioOS] Initialisation... ")
        # Initialiser le simulateur
        self.simulator = Simulator()
        self.status = "OFFLINE"
        self.config = config
        self.bsim_program = None
        
        # Configuration de l'interface matérielle
        self.hardware_enabled = config.get('hardware', {}).get('enabled', False)
        self.hardware_api_url = config.get('hardware', {}).get('api_url', 'http://localhost:5000')

    def compile_and_load_bsim(self, arbol_file_path: str):
        """Compile un fichier Arbol (.arbol) et charge le programme .bsim résultant."""
        print(f"[BioOS] Compilation du fichier Arbol: {arbol_file_path}...")
        
        if not os.path.exists(arbol_file_path):
            print(f"[BioOS ERROR] Fichier Arbol non trouvé: {arbol_file_path}")
            return

        with open(arbol_file_path, 'r') as f:
            code = f.read()

        try:
            output_json_path = arbol_file_path.replace('.arbol', '.bsim.json')
            subprocess.run([
                sys.executable, '-m', 'arbol.compiler.compiler', arbol_file_path
            ], check=True)
            with open(output_json_path, 'r') as jf:
                self.bsim_program = json.load(jf)
            print("[BioOS] Compilation Arbol réussie. Programme bsim chargé.")
        except Exception as e:
            print(f"[BioOS ERROR] Échec de la compilation Arbol: {e}")
            self.bsim_program = None

    def start(self):
        """Démarre le BioOS et initialise le simulateur."""
        print(f"[BioOS] Démarrage de la PQPE...")
        self.status = "BOOTING"
        
        if self.hardware_enabled:
            print(f"[BioOS] Vérification de la connexion matérielle ({self.hardware_api_url})...")
            try:
                response = requests.get(f"{self.hardware_api_url}/api/status", timeout=2)
                if response.status_code == 200:
                    print("[BioOS] Connexion matérielle établie.")
                else:
                    print(f"[BioOS WARNING] Matériel non prêt: {response.json()}")
            except Exception as e:
                print(f"[BioOS WARNING] Impossible de se connecter au Jetson Client: {e}")
                self.hardware_enabled = False

        print("[BioOS] Initialisation des modules biologiques et environnementaux...")
        self.status = "IDLE"
        print("[BioOS] Système PRÊT. En attente de programmes.")

    def _sync_hardware_instruction(self, command, params):
        """Synchronise une instruction logique avec le matériel réel."""
        if not self.hardware_enabled:
            return

        try:
            if command == 'STIMULUS_APPLY':
                stim = params.get('stimulus')
                args = params.get('arguments', {})
                if stim == 'light':
                    requests.post(f"{self.hardware_api_url}/api/actuators/light/set", json={
                        "channel_id": args.get('channel', 'DEFAULT'),
                        "intensity": args.get('intensity', 0.0),
                        "duration_ms": args.get('duration_ms', 1000)
                    }, timeout=1)
                elif stim == 'em_field':
                    requests.post(f"{self.hardware_api_url}/api/actuators/em_field/set", json={
                        "frequency_hz": args.get('frequency', 0.0),
                        "amplitude_v": args.get('amplitude', 0.0)
                    }, timeout=1)
            
            elif command == 'MEASURE':
                res = requests.get(f"{self.hardware_api_url}/api/sensors/electrode/read", 
                                 params={"channel_id": params.get('qubit')}, timeout=1)
                if res.status_code == 200:
                    print(f"[BioOS] Lecture matérielle reçue: {res.json().get('data')}")

        except Exception as e:
            print(f"[BioOS ERROR] Échec de synchronisation matérielle: {e}")

    def execute_bsim(self):
        """Exécute le programme bsim actuellement chargé."""
        if self.status != "IDLE":
            print("[BioOS ERROR] Le système n'est pas prêt. Veuillez démarrer l'OS.")
            return

        if self.bsim_program is None:
            print("[BioOS ERROR] Aucun programme bsim n'est chargé.")
            return

        self.status = "RUNNING"
        print("[BioOS] Exécution du programme bsim...")

        # Exécution et synchronisation matérielle
        for instruction in self.bsim_program['instructions']:
            self._sync_hardware_instruction(instruction.get('command'), instruction.get('params', {}))

        init_config = self.bsim_program['instructions'][0]['config']
        
        total_duration = 0
        for instruction in self.bsim_program['instructions']:
            if instruction['command'] == 'RUN_UNTIL':
                run_time = instruction['params']['time']
                if run_time > total_duration:
                    total_duration = run_time
        
        time_step = init_config.get('dt', 0.01)

        print(f"[BioOS] Durée totale de la simulation (bsim): {total_duration}s")

        # Exécuter directement le script BSIM via le simulateur unifié
        simulation_log = self.simulator.run_script(self.bsim_program['instructions'])

        # Sauvegarder le log de simulation
        results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '05_data', 'results'))
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        log_file_path = os.path.join(results_dir, 'simulation_log.json')
        with open(log_file_path, 'w') as f:
            json.dump(simulation_log, f, indent=2)
        print(f"[BioOS] Log de simulation sauvegardé dans: {log_file_path}")

        # Calculer et sauvegarder les métriques de convergence par run
        metrics = self._compute_bsim_metrics(self.bsim_program)
        metrics_path = os.path.join(results_dir, 'bsim_metrics.json')
        with open(metrics_path, 'w') as f:
            json.dump(metrics, f, indent=2)
        print(f"[BioOS] Métriques BSIM sauvegardées dans: {metrics_path}")

        # Générer une figure de convergence (fidélité/perte)
        try:
            import matplotlib.pyplot as plt
            runs = [m['run'] for m in metrics['trainer_runs']]
            fid = [m['fidelity'] for m in metrics['trainer_runs']]
            loss = [m['loss'] for m in metrics['trainer_runs']]

            plt.figure(figsize=(8, 5))
            plt.plot(runs, fid, marker='o', label='Fidélité (proxy)')
            plt.plot(runs, loss, marker='s', label='Perte (proxy)')
            plt.xlabel('Run (trainer)')
            plt.ylabel('Valeur')
            plt.title('Convergence PhytoQMML (depuis BSIM)')
            plt.grid(True, alpha=0.3)
            plt.legend()
            fig_path = os.path.join(results_dir, 'bsim_convergence.png')
            plt.tight_layout()
            plt.savefig(fig_path)
            print(f"[BioOS] Figure de convergence sauvegardée: {fig_path}")
        except Exception as e:
            print(f"[BioOS] Génération de figure ignorée ({e}).")

        print("[BioOS] Exécution du programme bsim terminée.")
        self.status = "IDLE"

    def stop(self):
        """Arrête la PQPE de manière contrôlée."""
        print("[BioOS] Arrêt de la PQPE...")
        self.status = "SHUTTING_DOWN"
        time.sleep(0.1)
        self.status = "OFFLINE"
        print("[BioOS] PQPE est maintenant hors ligne.")

    def _compute_bsim_metrics(self, program):
        """Calcule des métriques (proxy) de fidélité et perte par run trainer à partir des instructions BSIM."""
        instructions = program.get('instructions', [])
        init = instructions[0] if instructions else {}
        decoherence_rate = init.get('config', {}).get('quantum', {}).get('decoherence_rate', 0.05)

        # Regrouper les instructions détaillées par run 'phytoqmmml_trainer'
        groups = []
        current = None
        for instr in instructions[1:]:  # skip INITIALIZE
            cmd = instr.get('command')
            if cmd == 'RUN_START' and instr.get('params', {}).get('circuit') == 'phytoqmmml_trainer':
                if current:
                    groups.append(current)
                current = []
            elif current is not None and cmd in ['GATE_APPLY', 'MEASURE', 'RUN_UNTIL', 'STIMULUS_APPLY']:
                current.append(instr)
        if current:
            groups.append(current)

        trainer_metrics = []
        for idx, grp in enumerate(groups, start=1):
            ops = [i for i in grp if i.get('command') == 'GATE_APPLY']
            measures = [i for i in grp if i.get('command') == 'MEASURE']
            cnots = [i for i in ops if i.get('params', {}).get('gate') == 'cnot']

            complexity = len(ops) + 2 * len(cnots)
            # Fidélité proxy décroit avec dé-cohérence et complexité; légère amélioration par run
            import math
            fid_base = math.exp(-decoherence_rate * max(1, complexity))
            fidelity = min(1.0, fid_base * (1 + 0.02 * idx))
            loss = max(0.0, 1.0 - fidelity + 0.02 * len(measures))

            trainer_metrics.append({
                'run': idx,
                'complexity': complexity,
                'ops': len(ops),
                'measures': len(measures),
                'fidelity': fidelity,
                'loss': loss
            })

        # Compteurs globaux
        total_ops = sum(1 for i in instructions if i.get('command') == 'GATE_APPLY')
        total_measures = sum(1 for i in instructions if i.get('command') == 'MEASURE')
        total_runs = sum(1 for i in instructions if i.get('command') == 'RUN_START')

        return {
            'trainer_runs': trainer_metrics,
            'global_counts': {
                'total_ops': total_ops,
                'total_measures': total_measures,
                'total_runs': total_runs
            }
        }

if __name__ == '__main__':
    config_path = os.path.join(project_root, '03_unified_simulator', 'config', 'default.json')
    arbol_program_path = os.path.join(project_root, 'arbol', 'test_hardware.arbol')

    print(f"Chargement de la configuration depuis: {config_path}")
    print(f"Chargement du programme Arbol depuis: {arbol_program_path}")

    with open(config_path, 'r') as f:
        config = json.load(f)

    os_instance = BioOS(config=config)
    os_instance.start()
    os_instance.compile_and_load_bsim(arbol_program_path)
    os_instance.execute_bsim()
    os_instance.stop()
