import numpy as np

# Model parameters
params = {
    'k_prod_p700': 0.1,   # Production rate of P700
    'k_deg_p700': 0.02,   # Degradation rate of P700
    'K_light': 0.5,       # Michaelis-Menten constant for light activation
    'n_light': 2          # Hill coefficient for light activation
}

# Initial conditions
p_initial = [0.0] # Initial concentration of P700

def gene_regulation_model(p, t, light_intensity, params):
    """
    Defines the ODE for P700 gene regulation.
    p: array of concentrations [P700]
    t: time
    light_intensity: current light intensity (0 to 1)
    params: dictionary of model parameters
    """
    P700 = p[0]

    # Hill equation for light-activated production
    production_rate = params['k_prod_p700'] * (light_intensity**params['n_light']) / (params['K_light']**params['n_light'] + light_intensity**params['n_light'])

    # Degradation
    degradation_rate = params['k_deg_p700'] * P700

    # Rate of change
    dP700_dt = production_rate - degradation_rate

    return [dP700_dt]
