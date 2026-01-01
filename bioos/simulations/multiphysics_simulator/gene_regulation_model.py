
import numpy as np

# Paramètres du modèle (Alignés sur HAWRA_Full_Scientific_Paper.md Section 4.2)
params = {
    'K_light': 0.5,       # K_I : Constante de demi-saturation
    'V_max_synthesis': 0.2, # alpha_P : Taux de synthèse maximal
    'k_degradation': 0.05, # delta_P : Taux de dégradation
    'n_Hill': 2.0         # n : Coefficient de Hill (Cooperativité)
}

# Condition initiale
p_initial = [0.0] # Concentration initiale de P700

# Modèle d'équations différentielles (Cinétique de Hill)
def gene_regulation_model(p, t, light_intensity, params):
    P700 = p[0]
    
    # Équation pour la concentration de P700 avec cinétique de Hill
    # d[P]/dt = alpha_P * (I^n / (K^n + I^n)) - delta_P * [P]
    
    I_n = light_intensity ** params['n_Hill']
    K_n = params['K_light'] ** params['n_Hill']
    
    # Protection contre la division par zéro si I=0 et K=0 (peu probable ici)
    denom = K_n + I_n
    if denom == 0:
        term_production = 0
    else:
        term_production = I_n / denom
        
    synthesis_rate = params['V_max_synthesis'] * term_production
    degradation_rate = params['k_degradation'] * P700
    
    dP700_dt = synthesis_rate - degradation_rate
    
    return [dP700_dt]
