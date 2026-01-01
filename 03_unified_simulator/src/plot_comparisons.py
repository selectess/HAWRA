import json
import os
import matplotlib.pyplot as plt

def load_records(path):
    with open(path) as f:
        data = json.load(f)
    records = data['results'] if isinstance(data, dict) and 'results' in data else data
    return records

def plot_scenario(name, instr_path, sched_path, out_dir):
    instr = load_records(instr_path)
    sched = load_records(sched_path)

    # Align by index (assume uniform dt and comparable durations)
    n = min(len(instr), len(sched))
    t_instr = [r.get('time', i) for i, r in enumerate(instr[:n])]
    t_sched = [r.get('time', i) for i, r in enumerate(sched[:n])]

    li_instr = [r.get('light_intensity', 0.0) for r in instr[:n]]
    li_sched = [r.get('light_intensity', 0.0) for r in sched[:n]]
    p_instr = [r.get('p700_concentration', 0.0) for r in instr[:n]]
    p_sched = [r.get('p700_concentration', 0.0) for r in sched[:n]]
    c_instr = [r.get('coherence_factor', 0.0) for r in instr[:n]]
    c_sched = [r.get('coherence_factor', 0.0) for r in sched[:n]]

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(t_instr, li_instr, label='instructions')
    axes[0].plot(t_sched, li_sched, label='schedule', linestyle='--')
    axes[0].set_ylabel('light_intensity')
    axes[0].legend()

    axes[1].plot(t_instr, p_instr, label='instructions')
    axes[1].plot(t_sched, p_sched, label='schedule', linestyle='--')
    axes[1].set_ylabel('p700_concentration')
    axes[1].legend()

    axes[2].plot(t_instr, c_instr, label='instructions')
    axes[2].plot(t_sched, c_sched, label='schedule', linestyle='--')
    axes[2].set_ylabel('coherence_factor')
    axes[2].set_xlabel('time')
    axes[2].legend()

    plt.suptitle(f'Comparaison {name}: instructions vs schedule')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f'comparison_{name}.png')
    plt.tight_layout(rect=(0.0, 0.03, 1.0, 0.95))
    plt.savefig(out_path, dpi=150)
    print(out_path)

def plot_phytoqmmml(phyto_path, out_dir):
    with open(phyto_path) as f:
        data = json.load(f)
    qh = data.get('quantum_history', [])
    if not qh:
        print('No quantum_history found')
        return
    by_time = {}
    for ev in qh:
        t = ev.get('time', 0.0)
        k = 'op' if 'op' in ev else ('measure' if 'measure' in ev else ('run' if 'run' in ev else 'other'))
        by_time.setdefault(t, {'QUANTUM_OP': 0, 'MEASURE': 0, 'run': 0, 'other': 0})
        if k == 'op':
            by_time[t]['QUANTUM_OP'] += 1
        elif k == 'measure':
            by_time[t]['MEASURE'] += 1
        elif k == 'run':
            by_time[t]['run'] += 1
        else:
            by_time[t]['other'] += 1
    times = sorted(by_time.keys())
    ops = [by_time[t]['QUANTUM_OP'] for t in times]
    meas = [by_time[t]['MEASURE'] for t in times]
    runs = [by_time[t]['run'] for t in times]
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(times, ops, label='QUANTUM_OP')
    axes[0].plot(times, meas, label='MEASURE')
    axes[0].plot(times, runs, label='run')
    axes[0].set_ylabel('count')
    axes[0].set_title('Comptage par phase')
    axes[0].legend()

    # Observables Z par run (moyenne dernière valeur par timestamp)
    z_by_time = {}
    for ev in qh:
        if 'z_exp' in ev and isinstance(ev['z_exp'], list):
            t = ev.get('time', 0.0)
            series = [float(x) for x in ev['z_exp']]
            if series:
                z_by_time.setdefault(t, []).append(series[-1])
    z_times = sorted(z_by_time.keys())
    z_vals = [sum(z_by_time[t])/len(z_by_time[t]) if z_by_time[t] else 0.0 for t in z_times]
    axes[1].plot(z_times, z_vals, label='⟨Z⟩ par run', color='tab:orange')
    axes[1].set_xlabel('time')
    axes[1].set_ylabel('⟨Z⟩')
    axes[1].legend()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'phytoqmmml_quantum_counts_observables.png')
    plt.tight_layout(rect=(0.0, 0.03, 1.0, 0.95))
    plt.savefig(out_path, dpi=150)
    print(out_path)

def plot_spectral_response_demo(out_dir):
    from hawra_simulator.simulator import Simulator
    instructions = [
        {
            'command': 'INITIALIZE',
            'config': {
                'dt': 1.0,
                'env': {},
                'bio': {
                    'genes': {
                        'gA': {'basal_rate': 0.02, 'initial_expression': 0.05, 'light_sensitivity': 0.5, 'peak_wavelength': 650, 'peak_sigma': 25},
                        'gB': {'basal_rate': 0.02, 'initial_expression': 0.05, 'light_sensitivity': 0.5, 'peak_wavelength': 450, 'peak_sigma': 25}
                    },
                    'grn': {}
                }
            }
        },
        {
            'command': 'STIMULUS_APPLY',
            'params': {'stimulus': 'adaptive_light', 'arguments': {'intensity': 1.0, 'duration': 10, 'wavelength': 650}}
        },
        {'command': 'RUN_UNTIL', 'params': {'time': 12}}
    ]
    sim = Simulator({})
    res = sim.run_script(instructions)
    records = res['results']
    t = [r.get('time', i) for i, r in enumerate(records)]
    li = [r.get('light_intensity', 0.0) for r in records]
    gA = [r.get('genes', {}).get('gA', {}).get('expression_level', 0.0) for r in records]
    gB = [r.get('genes', {}).get('gB', {}).get('expression_level', 0.0) for r in records]
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    axes[0].plot(t, li, label='light_intensity (λ=650)')
    axes[0].set_ylabel('Light')
    axes[0].legend()
    axes[1].plot(t, gA, label='gA@650nm')
    axes[1].plot(t, gB, label='gB@450nm', linestyle='--')
    axes[1].set_ylabel('Expression')
    axes[1].set_xlabel('time')
    axes[1].legend()
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'gene_spectral_response.png')
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    print(out_path)

def plot_spectral_sweep(out_dir):
    from hawra_simulator.simulator import Simulator
    import csv
    wavelengths = list(range(450, 701, 25))
    rows = [('wavelength_nm', 'gA_expr', 'gB_expr')]
    gA_vals = []
    gB_vals = []
    for wl in wavelengths:
        instructions = [
            {
                'command': 'INITIALIZE',
                'config': {
                    'dt': 1.0,
                    'env': {},
                    'bio': {
                        'genes': {
                            'gA': {'basal_rate': 0.02, 'initial_expression': 0.05, 'light_sensitivity': 0.5, 'peak_wavelength': 650, 'peak_sigma': 25},
                            'gB': {'basal_rate': 0.02, 'initial_expression': 0.05, 'light_sensitivity': 0.5, 'peak_wavelength': 450, 'peak_sigma': 25}
                        },
                        'grn': {}
                    }
                }
            },
            {'command': 'STIMULUS_APPLY', 'params': {'stimulus': 'adaptive_light', 'arguments': {'intensity': 1.0, 'duration': 10, 'wavelength': wl}}},
            {'command': 'RUN_UNTIL', 'params': {'time': 12}}
        ]
        sim = Simulator({})
        res = sim.run_script(instructions)
        final_genes = res['final_state']['genes']
        gA = final_genes['gA']['expression_level']
        gB = final_genes['gB']['expression_level']
        rows.append((str(wl), f"{gA:.6f}", f"{gB:.6f}"))
        gA_vals.append(gA)
        gB_vals.append(gB)
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, 'spectral_sweep.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(wavelengths, gA_vals, label='gA (peak 650 nm)')
    ax.plot(wavelengths, gB_vals, label='gB (peak 450 nm)')
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Final expression')
    ax.legend()
    png_path = os.path.join(out_dir, 'spectral_sweep.png')
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    print(csv_path)
    print(png_path)

def main():
    base = os.path.dirname(os.path.dirname(__file__))
    out_dir = os.path.join(base, 'results')
    scenarios = [
        (
            'pulses_40_3_10',
            os.path.join(base, 'output_pulses_40_3_10.json'),
            os.path.join(base, 'output_from_arbol_pulses_40_3_10.json')
        ),
        (
            'pulses_40_1_5',
            os.path.join(base, 'output_pulses_40_1_5.json'),
            os.path.join(base, 'output_from_arbol_pulses_40_1_5.json')
        ),
        (
            'pulses_60_1_5',
            os.path.join(base, 'output_pulses_60_1_5.json'),
            os.path.join(base, 'output_from_arbol_pulses_60_1_5.json')
        ),
    ]
    for name, instr_path, sched_path in scenarios:
        if os.path.exists(instr_path) and os.path.exists(sched_path):
            plot_scenario(name, instr_path, sched_path, out_dir)
        else:
            print(f"Missing files for {name}: {instr_path} or {sched_path}")
    phyto = os.path.join(base, 'output_phytoqmmml_multi_runs.json')
    if os.path.exists(phyto):
        plot_phytoqmmml(phyto, out_dir)
    try:
        plot_spectral_response_demo(out_dir)
    except Exception as e:
        print(f"spectral demo failed: {e}")
    try:
        plot_spectral_sweep(out_dir)
    except Exception as e:
        print(f"spectral sweep failed: {e}")

if __name__ == '__main__':
    main()
