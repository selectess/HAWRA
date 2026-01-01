
import numpy as np

# Paramètres du modèle
params = {
    'K_light': 0.5,  # Constante de Michaelis-Menten pour la lumière
    'V_max_synthesis': 0.2, # Taux de synthèse maximal de P700
    'k_degradation': 0.05 # Taux de dégradation de P700
}

# Condition initiale
p_initial = [0.0] # Concentration initiale de P700

# Modèle d'équations différentielles
def gene_regulation_model(p, t, light_intensity, params):
    P700 = p[0]
    
    # Équation pour la concentration de P700
    synthesis_rate = params['V_max_synthesis'] * light_intensity / (params['K_light'] + light_intensity)
    degradation_rate = params['k_degradation'] * P700
    
    dP700_dt = synthesis_rate - degradation_rate
    
    return [dP700_dt]
