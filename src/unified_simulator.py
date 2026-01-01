
import os
import json
import numpy as np
from hawra_simulator.simulator import Simulator
from qutip import Qobj

def run_simulation_from_bsim(bsim_file_path):
    """
    Bridge function for the WebApp to run a simulation from a BSIM file.
    """
    simulator = Simulator(bsim_script=bsim_file_path)
    results = simulator.run()
    
    # Add a 'states' attribute to the results dictionary to satisfy app.py expectations
    # if it's not already there. We'll convert the final state to a Qobj.
    
    class SimulationResultWrapper:
        def __init__(self, data, simulator):
            self.data = data
            self.states = []
            
            # If we have a quantum state, convert it to a Qobj for Bloch sphere plotting
            if simulator.quantum_state:
                # Convert numpy array state vector to Qobj
                # Note: QuantumState.state_vector is a 1D numpy array
                q_state = Qobj(simulator.quantum_state.state_vector)
                self.states.append(q_state)
            
            # Also check quantum history for intermediate states if any were recorded
            # (Currently simulator.py only records ops, not states, but we could add it)
            
        def __getitem__(self, key):
            return self.data[key]

    return SimulationResultWrapper(results, simulator)
