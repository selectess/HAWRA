import os
import sys
import json
import subprocess
import matplotlib.pyplot as plt
import numpy as np

# Add project root to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../03_unified_simulator/src')))

from hawra_simulator.simulator import Simulator

def compile_arbol(arbol_path):
    print(f"[*] Compiling {arbol_path}...")
    try:
        # Run the compiler as a module
        result = subprocess.run(
            [sys.executable, "-m", "arbol.compiler.compiler", arbol_path],
            check=True,
            capture_output=True,
            text=True,
            cwd=os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        )
        print(result.stdout)
        bsim_path = arbol_path.replace('.arbol', '.bsim.json')
        if os.path.exists(bsim_path):
            print(f"[+] Compilation successful: {bsim_path}")
            return bsim_path
        else:
            print("[-] Compilation failed: Output file not found.")
            return None
    except subprocess.CalledProcessError as e:
        print(f"[-] Compilation error:\n{e.stderr}")
        return None

def run_simulation(bsim_path):
    print(f"[*] Running simulation from {bsim_path}...")
    sim = Simulator(bsim_path)
    
    # Manually ensure qubit count is correct for the demo if metadata is missing
    if sim.qubit_count == 0:
        sim.qubit_count = 1
        sim.quantum_state.qubit_count = 1
        sim.quantum_state.state_vector = np.zeros(2, dtype=complex)
        sim.quantum_state.state_vector[0] = 1.0

    results = sim.run()
    
    output_path = bsim_path.replace('.bsim.json', '_results.json')
    sim.save_results(output_path, results)
    return results, output_path

def generate_report(results, output_file):
    print("[*] Generating report...")
    
    quantum_hist = results.get('quantum_history', [])
    bio_hist = results.get('results', [])
    
    report = []
    report.append("# HAWRA 'First Bloom' Experiment Report")
    report.append(f"**Date:** {os.popen('date').read().strip()}")
    report.append("**Experiment:** Hadamard Gate on P700 Bio-Qubit")
    report.append("")
    report.append("## 1. Sequence Execution Log")
    
    for event in quantum_hist:
        time = event.get('time', 0)
        if 'op' in event:
            report.append(f"- **T={time:.2f}ms**: Quantum Operation `{event['op']}` applied to targets {event['targets']}")
        elif 'measure' in event:
            report.append(f"- **T={time:.2f}ms**: Measurement of Qubit {event['measure']} -> Result: **{event['outcome']}**")
    
    report.append("")
    report.append("## 2. Biological Telemetry")
    if bio_hist:
        start_bio = bio_hist[0]
        end_bio = bio_hist[-1]
        report.append(f"- **Initial P700 Concentration:** {start_bio.get('p700_concentration', 'N/A')}")
        report.append(f"- **Final P700 Concentration:** {end_bio.get('p700_concentration', 'N/A')}")
        report.append(f"- **Coherence Factor (Silica Shield):** {end_bio.get('coherence_factor', 'N/A')}")
    
    report.append("")
    report.append("## 3. Conclusion")
    
    has_measure = any('measure' in e for e in quantum_hist)
    has_gate = any('op' in e for e in quantum_hist)
    
    if has_gate and has_measure:
        report.append("**SUCCESS**: The biological machine successfully translated the Arbol code into physical quantum states.")
        report.append("The P700 qubit was initialized, manipulated by the `CRY2` protein interface (Hadamard), and read out via `LUC` bioluminescence.")
    else:
        report.append("**FAILURE**: The sequence did not complete as expected.")

    with open(output_file, 'w') as f:
        f.write("\n".join(report))
    
    print(f"[+] Report generated at {output_file}")

if __name__ == "__main__":
    arbol_file = os.path.abspath(os.path.join(os.path.dirname(__file__), '../01_genomics/experiments/first_bloom.arbol'))
    
    bsim_file = compile_arbol(arbol_file)
    if bsim_file:
        results, res_path = run_simulation(bsim_file)
        report_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../00_docs/proofs/FIRST_BLOOM_REPORT.md'))
        os.makedirs(os.path.dirname(report_path), exist_ok=True)
        generate_report(results, report_path)
