import numpy as np
from scipy.integrate import odeint
from .gene_regulation_model import gene_regulation_model, p_initial, params

class BiologicalEngine:
    def __init__(self, config):
        print("Initializing Biological Engine...")
        self.config = config
        self.p700_concentration = config.get('p700_initial', 0.0)
        
        # Override default params with config if provided
        self.params = {
            'K_light': config.get('K_light', 0.5),
            'V_max_synthesis': config.get('p700_synthesis_rate', 0.2),
            'k_degradation': config.get('p700_degradation_rate', 0.05),
            'n_Hill': config.get('n_Hill', 2.0)
        }

    def update(self, time, dt, env_state):
        light_intensity = env_state.get('light_intensity', 0)
        
        # Create a time array for the integration step
        t = [time, time + dt]
        
        # Solve the ODE for the current time step
        # The model function expects args in the order (light_intensity, params)
        solution = odeint(gene_regulation_model, [self.p700_concentration], t, args=(light_intensity, self.params))
        self.p700_concentration = solution[1][0]

        if time % 50 == 0 or time >= 490: # Reduce log verbosity but keep critical points
            print(f"Updating biology at t={time}. Light: {light_intensity}, P700: {self.p700_concentration:.4f}")

        return {
            'p700_concentration': self.p700_concentration
        }