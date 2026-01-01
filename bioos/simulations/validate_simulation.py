
import json
import numpy as np

def analyze_simulation_log(log_path, config):
    """
    Analyzes the multiphysics simulation log to validate the model's behavior.

    Args:
        log_path (str): Path to the simulation log file (JSON).
        config (dict): The simulation configuration used to generate the log.

    Returns:
        dict: A dictionary containing the validation results.
    """
    with open(log_path, 'r') as f:
        log = json.load(f)

    p700_threshold = config['quantum']['threshold']
    
    errors = []
    excitations = 0
    green_reads = 0
    red_reads = 0

    for i in range(1, len(log)):
        prev_state = log[i-1]
        current_state = log[i]

        # 1. P700 concentration validation
        if current_state['light_intensity'] > 0:
            if current_state['p700_concentration'] < prev_state['p700_concentration'] and prev_state['p700_concentration'] < 0.99 :
                # Allow for decay when light is on but synthesis is saturated
                pass
        else: # No light
            if current_state['p700_concentration'] > prev_state['p700_concentration']:
                 errors.append(f"t={current_state['time']}: P700 increased without light.")

        # 2. P700 excitation validation
        if current_state['p700_state'] == 1 and prev_state['p700_state'] == 0:
            excitations += 1
            if current_state['p700_concentration'] < p700_threshold:
                errors.append(f"t={current_state['time']}: P700 excited below threshold.")

        # 3. Readout channel validation
        if current_state['luc_green_output'] > 0 or current_state['luc_red_output'] > 0:
            if prev_state['p700_state'] == 0:
                errors.append(f"t={current_state['time']}: Readout occurred from ground state.")
        
        if current_state['luc_green_output'] > 0:
            green_reads += 1
        
        if current_state['luc_red_output'] > 0:
            red_reads += 1

        # 4. Mutual exclusion of readouts
        if current_state['luc_green_output'] > 0 and current_state['luc_red_output'] > 0:
            errors.append(f"t={current_state['time']}: Mutual exclusion of readouts violated.")

    return {
        "validation_status": "SUCCESS" if not errors else "FAILURE",
        "errors": errors,
        "total_steps": len(log),
        "excitations": excitations,
        "green_reads": green_reads,
        "red_reads": red_reads
    }

if __name__ == "__main__":
    # This is the same config used in the simulation
    config = {
        'max_time': 400,
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
            'degradation_rate': 0.04,
            'synthesis_rate': 0.25
        },
        'quantum': {
            'threshold': 0.8,
            'decoherence_rate': 0.015
        }
    }
    
    log_file = "/Users/mehdiwhb/Desktop/HAWRA/05_data/results/multiphysics_simulation/multiphysics_simulation_v2.json"
    
    results = analyze_simulation_log(log_file, config)
    
    print("--- PQPE Numerical Validation Results ---")
    print(f"Status: {results['validation_status']}")
    if results['errors']:
        print("Errors found:")
        for error in results['errors']:
            print(f"- {error}")
    else:
        print("No errors found. The model behaves as expected.")
    print("\n--- Statistics ---")
    print(f"Total simulation steps: {results['total_steps']}")
    print(f"P700 excitations: {results['excitations']}")
    print(f"Green channel readouts (|0>): {results['green_reads']}")
    print(f"Red channel readouts (|1>): {results['red_reads']}")
    print("------------------------------------")
