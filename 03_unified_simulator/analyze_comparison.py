
import json
import numpy as np

def calculate_final_state_probabilities(final_state):
    return [s['real']**2 + s['imag']**2 for s in final_state]

def generate_report(clean_data, noisy_data):
    report = """
# Rapport de Comparaison de Simulation

Ce rapport compare les résultats de deux simulations : une simulation "clean" (sans bruit) et une simulation "noisy" (avec bruit).

## 1. Paramètres de Simulation

| Paramètre | Simulation Clean | Simulation Noisy |
|---|---|---|
| Température | {clean_temp:.2f}°C | {noisy_temp:.2f}°C |
| Intensité Lumineuse (max) | {clean_light:.2f} | {noisy_light:.2f} |

## 2. Métriques de Performance Biologique

| Métrique | Simulation Clean | Simulation Noisy |
|---|---|---|
| Concentration P700 (finale) | {clean_p700:.4f} | {noisy_p700:.4f} |
| Stabilité (moyenne) | {clean_stability:.4f} | {noisy_stability:.4f} |
| Cohérence (moyenne) | {clean_coherence:.4f} | {noisy_coherence:.4f} |

## 3. Analyse de l'État Final

Les distributions de probabilité de l'état final pour les deux simulations sont présentées ci-dessous.

"""

    # Extract parameters
    clean_temp = clean_data['environment_history'][0]['temperature']
    noisy_temp = noisy_data['environment_history'][0]['temperature']
    clean_light = max(e['light_intensity'] for e in clean_data['environment_history'])
    noisy_light = max(e['light_intensity'] for e in noisy_data['environment_history'])

    # Extract biological metrics
    clean_p700 = clean_data['biological_system_history'][-1]['p700_concentration']
    noisy_p700 = noisy_data['biological_system_history'][-1]['p700_concentration']
    clean_stability = np.mean([b['stability'] for b in clean_data['biological_system_history']])
    noisy_stability = np.mean([b['stability'] for b in noisy_data['biological_system_history']])
    clean_coherence = np.mean([b['coherence_factor'] for b in clean_data['biological_system_history']])
    noisy_coherence = np.mean([b['coherence_factor'] for b in noisy_data['biological_system_history']])

    return report.format(
        clean_temp=clean_temp,
        noisy_temp=noisy_temp,
        clean_light=clean_light,
        noisy_light=noisy_light,
        clean_p700=clean_p700,
        noisy_p700=noisy_p700,
        clean_stability=clean_stability,
        noisy_stability=noisy_stability,
        clean_coherence=clean_coherence,
        noisy_coherence=noisy_coherence
    )

def main():
    with open('output_clean.json', 'r') as f:
        clean_data = json.load(f)

    with open('output_noisy.json', 'r') as f:
        noisy_data = json.load(f)

    report = generate_report(clean_data, noisy_data)
    print(report)

    # You can also save the report to a file
    with open('comparison_report.md', 'w') as f:
        f.write(report)

if __name__ == '__main__':
    main()
