import json
from .biological_engine import BiologicalEngine
from .quantum_engine import QuantumEngine
from .environment_engine import EnvironmentEngine

class MultiphysicsSimulator:
    def __init__(self, config):
        print("Initializing Multiphysics Simulator...")
        self.config = config
        self.bio_engine = BiologicalEngine(config.get('bio', {}))
        self.quantum_engine = QuantumEngine(config.get('quantum', {}))
        self.env_engine = EnvironmentEngine(config.get('env', {}))
        self.time = 0
        self.max_time = config.get('max_time', 100)
        self.dt = config.get('dt', 1)
        self.log = []

    def run(self):
        print("Starting simulation loop.")
        while self.time < self.max_time:
            self.step()
        print("Simulation finished.")
        return self.log

    def run_until(self, time):
        """Exécute la simulation jusqu'à un temps donné."""
        while self.time < time:
            self.step()

    def run_and_plot_steps(self, time, frame_dir):
        """Exécute la simulation et sauvegarde un graphique à chaque étape."""
        frame_number = 0
        while self.time < time:
            self.step()
            frame_path = f"{frame_dir}/frame_{frame_number:04d}.png"
            self.plot_step(frame_path)
            frame_number += 1

    def plot_step(self, output_path):
        """Génère un graphique de l'état actuel de la simulation."""
        try:
            import matplotlib.pyplot as plt
            times = [s['time'] for s in self.log]
            light = [s['light_intensity'] for s in self.log]
            p700 = [s['p700_concentration'] for s in self.log]
            luc_green = [s['luc_green_output'] for s in self.log]
            luc_red = [s['luc_red_output'] for s in self.log]

            fig, ax1 = plt.subplots(figsize=(12, 6))

            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('P700 Concentration', color='tab:blue')
            ax1.plot(times, p700, color='tab:blue', label='P700 Concentration')
            ax1.tick_params(axis='y', labelcolor='tab:blue')
            ax1.set_ylim(0, 1) # Fix Y-axis for consistency

            ax2 = ax1.twinx()
            ax2.set_ylabel('Signal')
            ax2.plot(times, light, color='orange', linestyle='--', label='Light Stimulus')
            ax2.bar(times, luc_green, width=self.dt, color='g', alpha=0.6, label='Canal Vert |0> (Stable)')
            ax2.bar(times, luc_red, width=self.dt, color='r', alpha=0.6, label='Canal Rouge |1> (Instable)')
            ax2.tick_params(axis='y')
            ax2.set_ylim(0, 1.1) # Fix Y-axis for consistency

            fig.tight_layout()
            plt.title(f'HAWRA Simulation - Time: {self.time:.2f}s')
            plt.legend()
            plt.savefig(output_path)
            plt.close(fig) # Close the figure to free memory

        except ImportError:
            print("\nMatplotlib not found. Skipping plot.")



    def step(self):
        # 1. Get light intensity from environment
        env_state = self.env_engine.update(self.time, self.dt)
        light_intensity = env_state['light_intensity']

        # 2. Update biological state
        self.bio_engine.update(self.time, self.dt, env_state)
        
        # 3. Update quantum state
        self.quantum_engine.update_state(self.bio_engine.p700_concentration)
        quantum_state = self.quantum_engine.get_state()

        # 4. Log current state
        self.log.append({
            'time': self.time,
            'light_intensity': light_intensity,
            'p700_concentration': self.bio_engine.p700_concentration,
            'p700_state': quantum_state['p700_state'],
            'luc_green_output': quantum_state['luc_green_output'],
            'luc_red_output': quantum_state['luc_red_output']
        })

        # 5. Advance time
        self.time += self.dt

    def plot_results(self, output_path):
        """Génère un graphique des résultats de la simulation."""
        try:
            import matplotlib.pyplot as plt
            times = [s['time'] for s in self.log]
            light = [s['light_intensity'] for s in self.log]
            p700 = [s['p700_concentration'] for s in self.log]
            luc_green = [s['luc_green_output'] for s in self.log]
            luc_red = [s['luc_red_output'] for s in self.log]

            fig, ax1 = plt.subplots(figsize=(12, 6))

            ax1.set_xlabel('Time (s)')
            ax1.set_ylabel('P700 Concentration', color='tab:blue')
            ax1.plot(times, p700, color='tab:blue', label='P700 Concentration')
            ax1.tick_params(axis='y', labelcolor='tab:blue')

            ax2 = ax1.twinx()
            ax2.set_ylabel('Signal')
            ax2.plot(times, light, color='orange', linestyle='--', label='Light Stimulus')
            # Utiliser des barres pour les événements de lecture discrets
            ax2.bar(times, luc_green, width=self.dt, color='g', alpha=0.6, label='Canal Vert |0> (Stable)')
            ax2.bar(times, luc_red, width=self.dt, color='r', alpha=0.6, label='Canal Rouge |1> (Instable)')
            ax2.tick_params(axis='y')

            fig.tight_layout()
            plt.title('HAWRA Genetic-Compliant Simulation')
            plt.legend()
            plt.savefig(output_path)
            print(f"\nExecution plot saved to {output_path}")

        except ImportError:
            print("\nMatplotlib not found. Skipping plot.")

if __name__ == "__main__":
    # Configuration de la simulation
    config = {
        'max_time': 400, # Temps de simulation plus long
        'dt': 0.5,
        'env': {
            'pulse_configs': [
                {'start': 10, 'end': 20, 'intensity': 1.0},
                {'start': 50, 'end': 55, 'intensity': 0.8},
                {'start': 90, 'end': 92, 'intensity': 1.0},
                {'start': 120, 'end': 140, 'intensity': 0.6},
                {'start': 160, 'end': 170, 'intensity': 1.0},
                {'start': 200, 'end': 220, 'intensity': 0.9},
                {'start': 250, 'end': 270, 'intensity': 0.7},
                {'start': 300, 'end': 310, 'intensity': 1.0},
                {'start': 340, 'end': 360, 'intensity': 0.5}
            ]
        },
        'bio': {
            'p700_initial': 0.0,
            'degradation_rate': 0.04, # Dégradation légèrement plus lente
            'synthesis_rate': 0.25 # Synthèse légèrement plus rapide
        },
        'quantum': {
            'threshold': 0.8, # Seuil légèrement plus bas
            'decoherence_rate': 0.015 # Décohérence plus faible
        }
    }

    # Initialiser et exécuter le simulateur
    simulator = MultiphysicsSimulator(config)
    log = simulator.run()

    # Sauvegarder les résultats
    output_dir = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/multiphysics_simulation"
    import os
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Sauvegarder le log
    log_path = os.path.join(output_dir, "multiphysics_simulation_v2.json")
    with open(log_path, 'w') as f:
        json.dump(log, f, indent=2)
    print(f"Simulation log saved to {log_path}")

    # Générer le graphique final
    plot_path = os.path.join(output_dir, "multiphysics_simulation_v2.png")
    simulator.plot_results(plot_path)