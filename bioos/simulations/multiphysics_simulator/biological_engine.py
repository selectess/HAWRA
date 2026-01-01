import numpy as np
from scipy.integrate import odeint
from .gene_regulation_model import gene_regulation_model, p_initial, params

class BiologicalEngine:
    def __init__(self, config):
        print("Initializing Biological Engine...")
        self.config = config
        self.p700_concentration = p_initial[0]

    def update(self, time, dt, env_state):
        light_intensity = env_state.get('light_intensity', 0)
        
        # Create a time array for the integration step
        t = [time, time + dt]
        
        # Solve the ODE for the current time step
        # The model function expects args in the order (light_intensity, params)
        solution = odeint(gene_regulation_model, [self.p700_concentration], t, args=(light_intensity, params))
        self.p700_concentration = solution[1][0]

        print(f"Updating biology at t={time}. P700 concentration: {self.p700_concentration:.4f}")

        return {
            'p700_concentration': self.p700_concentration
        }