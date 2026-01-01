import json
import os
import matplotlib.pyplot as plt
import numpy as np

def load_simulation_output(path):
    with open(path, 'r') as f:
        return json.load(f)

def plot_final_state_probabilities(clean_data, noisy_data, out_path):
    # Process clean data
    final_state_complex_clean = [d['real'] + 1j * d['imag'] for d in clean_data['final_state']]
    probabilities_clean = np.abs(np.array(final_state_complex_clean))**2
    
    # Process noisy data
    final_state_complex_noisy = [d['real'] + 1j * d['imag'] for d in noisy_data['final_state']]
    probabilities_noisy = np.abs(np.array(final_state_complex_noisy))**2
    
    num_qubits = int(np.log2(len(probabilities_clean)))
    state_labels = [format(i, f'0{num_qubits}b') for i in range(len(probabilities_clean))]

    x = np.arange(len(state_labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(15, 7))
    rects1 = ax.bar(x - width/2, probabilities_clean, width, label='Sans Bruit (35°C)', color='#4e79a7')
    rects2 = ax.bar(x + width/2, probabilities_noisy, width, label='Avec Bruit (45°C)', color='#f28e2b')

    ax.set_title("Comparaison des Probabilités de l'état final (Avec vs Sans Bruit)")
    ax.set_xlabel('État de la base de calcul')
    ax.set_ylabel('Probabilité')
    ax.set_xticks(x)
    ax.set_xticklabels(state_labels, rotation=90)
    ax.legend()

    fig.tight_layout()
    plt.savefig(out_path)
    print(f"Figure sauvegardée: {out_path}")

if __name__ == '__main__':
    base = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    clean_output_path = os.path.join(base, '03_unified_simulator', 'output_clean.json')
    noisy_output_path = os.path.join(base, '03_unified_simulator', 'output_noisy.json')
    out_path = os.path.join(base, '03_unified_simulator', 'results', 'final_state_comparison.png')
    
    clean_data = load_simulation_output(clean_output_path)
    noisy_data = load_simulation_output(noisy_output_path)
    
    plot_final_state_probabilities(clean_data, noisy_data, out_path)
