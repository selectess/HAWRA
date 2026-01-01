import os
import json
import matplotlib.pyplot as plt
import numpy as np

# Essayer d'abord la visualisation quantique depuis validation_results.json,
# sinon visualiser la dynamique GRN depuis simulation_results.json
base = '/Users/mehdiwhb/Desktop/HAWRA/05_data/results'
validation_path = os.path.join(base, 'validation_results.json')
simulation_path = os.path.join(base, 'simulation_results.json')

if os.path.exists(validation_path):
    from qutip import Bloch, Qobj
    from matplotlib.animation import FuncAnimation
    with open(validation_path, 'r') as f:
        results = json.load(f)
    quantum_history = results.get('quantum_history', [])
    if quantum_history:
        states_raw = [entry['state_vector'] for entry in quantum_history]
        states = []
        for s in states_raw:
            complex_vector = np.array([c['real'] + 1j * c['imag'] for c in s])
            states.append(complex_vector)
        fig = plt.figure()
        b = Bloch()
        def animate(i):
            b.clear()
            state_vector = states[i]
            if len(state_vector) >= 2:
                q_state = Qobj(state_vector[[0, 1]])
                b.add_states(q_state.unit())
            b.make_sphere()
            return b.fig
        anim = FuncAnimation(fig, animate, frames=len(states), repeat=False)
        anim.save(os.path.join(base, 'bloch_animation.gif'), writer='imagemagick', fps=5)
        print('Animation de la sphère de Bloch sauvegardée.')
    else:
        print('validation_results.json présent mais sans historique quantique, bascule vers GRN.')

elif os.path.exists(simulation_path):
    with open(simulation_path, 'r') as f:
        data = json.load(f)
    bio_hist = data.get('results', [])
    env_hist = data.get('environment_history', [])
    times = [entry.get('time', None) for entry in env_hist]
    # Accumuler les courbes par gène
    gene_series = {}
    for snapshot in bio_hist:
        genes = snapshot.get('genes', {})
        for gid, gvals in genes.items():
            gene_series.setdefault(gid, []).append(gvals.get('expression_level', 0.0))
    # Tracer
    fig, ax = plt.subplots(figsize=(8, 5))
    for gid, series in gene_series.items():
        if times and len(times) == len(series):
            ax.plot(times, series, label=gid)
        else:
            ax.plot(series, label=gid)
    ax.set_xlabel('Temps')
    ax.set_ylabel('Expression (normalisée)')
    ax.set_title('Dynamique du réseau de régulation génique')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left')
    ax2 = ax.twinx()
    li = [e.get('light_intensity', 0.0) for e in env_hist]
    if times and li and len(times) == len(li):
        ax2.plot(times, li, color='black', linestyle='--', alpha=0.6, label='Lumière')
        ax2.set_ylabel('Intensité lumineuse')
        ax2.legend(loc='upper right')
        # Zones ombrées pour les pulses
        in_pulse = False
        start_t = None
        labeled = False
        for i in range(1, len(times)):
            prev = li[i-1]
            curr = li[i]
            if not in_pulse and prev == 0 and curr > 0:
                in_pulse = True
                start_t = times[i]
            elif in_pulse and prev > 0 and curr == 0:
                end_t = times[i]
                if labeled:
                    ax.axvspan(start_t, end_t, color='yellow', alpha=0.2)
                else:
                    ax.axvspan(start_t, end_t, color='yellow', alpha=0.2, label='Pulse')
                    labeled = True
                in_pulse = False
                start_t = None
        if in_pulse and start_t is not None:
            end_t = times[-1]
            if labeled:
                ax.axvspan(start_t, end_t, color='yellow', alpha=0.2)
            else:
                ax.axvspan(start_t, end_t, color='yellow', alpha=0.2, label='Pulse')
                labeled = True
    fig.tight_layout()
    # Exporter CSV des statistiques par gène
    import csv
    stats_path = os.path.join(base, 'grn_stats.csv')
    with open(stats_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['gene', 'max', 'min', 'peak_time'])
        for gid, series in gene_series.items():
            if not series:
                writer.writerow([gid, '', '', ''])
                continue
            vmax = max(series)
            vmin = min(series)
            idx = series.index(vmax)
            peak_time = times[idx] if times and len(times) == len(series) else idx
            writer.writerow([gid, vmax, vmin, peak_time])

    out_png = os.path.join(base, 'grn_dynamics.png')
    fig.savefig(out_png)
    print(f'Dynamique GRN sauvegardée: {out_png}')
    print(f'Statistiques GRN sauvegardées: {stats_path}')
else:
    print('Aucun fichier de résultats trouvé pour la visualisation.')
