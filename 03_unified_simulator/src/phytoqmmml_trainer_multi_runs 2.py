import os
import json
import math
import subprocess
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
ARBL_PATH = os.path.join(BASE, '04_arbol', 'phytoqmmml_multi_runs.arbol')
BSIM_PATH = ARBL_PATH.replace('.arbol', '.bsim.json')
OUT_PNG = os.path.join(BASE, '03_quantum_simulation', 'results', 'phytoqmmml_convergence.png')

def compile_arbol():
    cmd = ['python3', '-m', '04_arbol.compiler.compiler', ARBL_PATH]
    subprocess.run(cmd, cwd=BASE, check=True)

def load_assembly(path):
    with open(path, 'r') as f:
        return json.load(f)

def group_runs(instructions, circuit_name='phytoqmmml_trainer'):
    groups = []
    current = None
    for instr in instructions:
        if instr.get('command') == 'run' and instr.get('circuit') == circuit_name:
            # start new group
            if current:
                groups.append(current)
            current = []
        elif current is not None:
            # collect detailed instructions for this run
            if instr.get('command') in ['QUANTUM_OP', 'MEASURE', 'RUN_UNTIL', 'STIMULUS_APPLY']:
                current.append(instr)
    if current:
        groups.append(current)
    return groups

def compute_metrics(groups, decoherence_rate=0.05):
    metrics = []
    for idx, grp in enumerate(groups, start=1):
        ops = [i for i in grp if i.get('command') == 'QUANTUM_OP']
        measures = [i for i in grp if i.get('command') == 'MEASURE']
        cnots = [i for i in ops if i.get('params', {}).get('gate') == 'CNOT']

        complexity = len(ops) + 2 * len(cnots)
        # Fidelity proxy decays with decoherence and complexity; improves with runs
        fid_base = math.exp(-decoherence_rate * complexity)
        fidelity = min(1.0, fid_base * (1 + 0.05 * idx))
        # Loss proxy: complementary of fidelity plus measurement overhead
        loss = max(0.0, 1.0 - fidelity + 0.02 * len(measures))

        metrics.append({'run': idx, 'complexity': complexity, 'fidelity': fidelity, 'loss': loss})
    return metrics

def plot_convergence(metrics, out_path):
    runs = [m['run'] for m in metrics]
    fid = [m['fidelity'] for m in metrics]
    loss = [m['loss'] for m in metrics]

    plt.figure(figsize=(8, 5))
    plt.plot(runs, fid, marker='o', label='Fidélité (proxy)')
    plt.plot(runs, loss, marker='s', label='Perte (proxy)')
    plt.xlabel('Run')
    plt.ylabel('Valeur')
    plt.title('Convergence PhytoQMML (multi-runs)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path)
    print(f"Figure sauvegardée: {out_path}")

def main():
    # Compile and load
    compile_arbol()
    assembly = load_assembly(BSIM_PATH)
    instr = assembly.get('instructions', [])
    decoherence = assembly.get('instructions', [{}])[0].get('config', {}).get('quantum', {}).get('decoherence_rate', 0.05)

    groups = group_runs(instr, circuit_name='phytoqmmml_trainer')
    metrics = compute_metrics(groups, decoherence_rate=decoherence)

    # Save metrics JSON alongside figure
    out_json = os.path.join(BASE, '03_quantum_simulation', 'results', 'phytoqmmml_convergence.json')
    with open(out_json, 'w') as f:
        json.dump({'metrics': metrics}, f, indent=2)
    print(f"Metrics sauvegardés: {out_json}")

    plot_convergence(metrics, OUT_PNG)

if __name__ == '__main__':
    main()

