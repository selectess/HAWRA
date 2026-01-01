
import os
import imageio
import json
from ..simulations.multiphysics_simulator.simulator import MultiphysicsSimulator

class HAWRACore:
    def __init__(self):
        self.simulator = None
        self.program = None
        self.frame_dir = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/temp_gif_frames"
        if not os.path.exists(self.frame_dir):
            os.makedirs(self.frame_dir)

    def load_program_from_file(self, file_path):
        """Charge un programme d'assembly quantique depuis un fichier JSON."""
        with open(file_path, 'r') as f:
            self.program = json.load(f)
        print(f"Programme chargé depuis {file_path}")

    def execute_program_and_generate_frames(self):
        """Exécute le programme chargé et génère les images pour le GIF."""
        if not self.program:
            print("Erreur: Aucun programme chargé.")
            return

        for instruction in self.program['instructions']:
            command = instruction.get('command')
            if command == "INITIALIZE":
                self.simulator = MultiphysicsSimulator(instruction.get('config'))
                print("Simulateur initialisé.")
            elif command == "RUN_UNTIL":
                if self.simulator:
                    time = instruction.get('params', {}).get('time')
                    self.simulator.run_and_plot_steps(time, self.frame_dir)
                else:
                    print("Erreur: Le simulateur n'est pas initialisé.")
            else:
                print(f"Commande inconnue: {command}")

    def create_gif(self, output_path, fps=10):
        """Crée un GIF à partir des images générées."""
        images = []
        for i in sorted(os.listdir(self.frame_dir)):
            if i.endswith('.png'):
                file_path = os.path.join(self.frame_dir, i)
                images.append(imageio.imread(file_path))
        
        imageio.mimsave(output_path, images, fps=fps)
        print(f"GIF sauvegardé dans {output_path}")

    def cleanup_frames(self):
        """Supprime les images temporaires."""
        for file_name in os.listdir(self.frame_dir):
            file_path = os.path.join(self.frame_dir, file_name)
            os.remove(file_path)
        os.rmdir(self.frame_dir)
        print("Images temporaires supprimées.")

if __name__ == '__main__':
    core = HAWRACore()
    
    core.load_program_from_file('00_docs/arbol_language/example.qasm.json')
    
    core.execute_program_and_generate_frames()
    
    output_gif = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/hawra_poc_execution_v2.gif"
    core.create_gif(output_gif, fps=20)
    
    core.cleanup_frames()
