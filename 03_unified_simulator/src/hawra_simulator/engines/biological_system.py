import numpy as np

class BiologicalSystem:
    """
    Represents the biological component of the HAWRA simulator.

    This class is responsible for managing the state of the biological system,
    including the gene regulatory network (GRN), the P700 concentration, and
    the overall stability of the system.
    """
    def __init__(self, config):
        """
        Initializes the biological system.

        Args:
            config (dict): A dictionary containing the configuration for the
                biological system.
        """
        self.config = config
        self.p700_concentration = config.get('initial_p700', 1.0)
        self.stability = 1.0
        self.coherence_factor = 1.0
        
        # Initialiser les gènes et leurs niveaux d'expression
        self.genes = {gene['id']: {'expression_level': gene.get('initial_expression', 0.0)} 
                      for gene in config.get('genes', [])}
        
        # Stocker les paramètres des gènes pour un accès facile
        self.gene_params = {gene['id']: gene for gene in config.get('genes', [])}
        
        # Réseau de régulation génique (GRN)
        self.grn = config.get('grn', {})
        
        self.synthesis_rate = config.get('synthesis_rate', 0.05)

    def update_state(self, temperature, dt=1.0):
        optimal = float(self.config.get('optimal_temp', 35.0))
        sigma = float(self.config.get('temp_sigma', 2.0))
        delta = float(temperature) - optimal
        self.stability = float(np.exp(- (delta ** 2) / (2 * sigma ** 2)))
        self.stability = float(np.clip(self.stability, 0.0, 1.0))
        self.coherence_factor = self.stability
        return self.get_state()

    def apply_stimulus(self, stimulus_params):
        """
        Applies a stimulus to the biological system.

        Args:
            stimulus_params (dict): A dictionary containing the parameters for
                the stimulus.
        """
        stimulus_type = stimulus_params.get('stimulus')
        target = stimulus_params.get('target')
        args = stimulus_params.get('arguments', {})

        if stimulus_type == 'inhibitor' and target in self.gene_params:
            concentration = args.get('concentration', 0.1)
            self.gene_params[target]['basal_rate'] *= (1 - concentration)
            print(f"[Stimulus] Applied inhibitor to {target}, new basal rate: {self.gene_params[target]['basal_rate']}")
        elif stimulus_type == 'chemical' and target in self.gene_params:
            effect = str(args.get('effect', 'downregulate')).lower()
            magnitude = float(args.get('magnitude', 0.1))
            if effect == 'upregulate':
                self.gene_params[target]['basal_rate'] *= (1 + magnitude)
            else:
                self.gene_params[target]['basal_rate'] *= max(0.0, (1 - magnitude))
            print(f"[Stimulus] Applied chemical to {target}, new basal rate: {self.gene_params[target]['basal_rate']}")

    def get_state(self):
        """
        Returns the current state of the biological system.

        Returns:
            dict: A dictionary containing the current state of the biological
                system.
        """
        # Inclure les niveaux d'expression des gènes dans l'état
        state = {
            'p700_concentration': self.p700_concentration,
            'stability': self.stability,
            'coherence_factor': self.coherence_factor,
            'genes': self.genes
        }
        return state

    def update(self, dt, light_intensity, light_wavelength=None):
        """
        Updates the state of the biological system for a given time step,
        incorporating polygenic interactions with Hill kinetics.

        Args:
            dt (float): The time step.
            light_intensity (float): The current light intensity.
        """
        new_expression_levels = {}

        for gene_id, params in self.gene_params.items():
            # Basal and environmental influence
            basal_rate = params.get('basal_rate', 0.01)
            light_sensitivity = params.get('light_sensitivity', 0)
            peak = params.get('peak_wavelength')
            sigma = float(params.get('peak_sigma', 30.0))
            if light_wavelength is not None and peak is not None:
                spectral_weight = float(np.exp(- ((float(light_wavelength) - float(peak)) ** 2) / (2 * sigma ** 2)))
            else:
                spectral_weight = 1.0
            effective_light = light_intensity * spectral_weight
            environmental_influence = basal_rate + light_sensitivity * effective_light

            # Regulatory influence from the GRN
            activation_influence = 0.0
            repression_influence = 1.0  # Multiplicative factor

            if gene_id in self.grn:
                regulators = self.grn[gene_id]

                # Handle old dictionary-based GRN format for backward compatibility
                if isinstance(regulators, dict):
                    temp_regulators = []
                    for reg_id, weight in regulators.items():
                        reg_type = "activator" if weight > 0 else "repressor"
                        # This simple conversion doesn't support Hill params for old format
                        temp_regulators.append({"id": reg_id, "type": reg_type, "weight": abs(weight)})
                    regulators = temp_regulators

                for reg_info in regulators:
                    regulator_id = reg_info["id"]
                    if regulator_id not in self.genes:
                        continue # Skip if regulator gene is not defined

                    regulator_level = self.genes[regulator_id]['expression_level']
                    
                    reg_type = reg_info.get("type", "activator")
                    weight = reg_info.get("weight", 1.0)
                    n = reg_info.get('hill_coefficient', 1)
                    k = reg_info.get('half_max_concentration', 0.5)

                    if reg_type == "activator":
                        # Hill function for activation, contribution is additive
                        activation = (weight * (regulator_level**n)) / (k**n + regulator_level**n)
                        activation_influence += activation
                    elif reg_type == "repressor":
                        # Hill function for repression, contribution is multiplicative
                        repression = (k**n) / (k**n + regulator_level**n)
                        # The final repression factor is a product of all repressor effects
                        repression_influence *= repression

            # Total synthesis rate combines basal, environmental, and regulatory influences
            synthesis_rate = (environmental_influence + activation_influence) * repression_influence

            # Update expression level using the synthesis-degradation model
            current_level = self.genes[gene_id]['expression_level']
            degradation_rate = params.get('degradation_rate', 0.05)
            delta_expression = (synthesis_rate - degradation_rate * current_level) * dt
            
            new_level = current_level + delta_expression
            new_expression_levels[gene_id] = {'expression_level': np.clip(new_level, 0, 1)}

        # Mettre à jour tous les niveaux d'expression des gènes en même temps
        self.genes = new_expression_levels

        # 2. Mettre à jour la stabilité et la concentration de P700
        # La stabilité est affectée par l'expression des gènes (par exemple, un gène de stabilité)
        stability_influence = 0.0
        if 'gene_stability' in self.genes: # Supposons qu'un gène affecte la stabilité
            stability_influence = self.genes['gene_stability']['expression_level'] * 0.1

        stability_decay = 0.01
        self.stability += (stability_influence - stability_decay) * dt
        self.stability = np.clip(self.stability, 0, 1)

        # La concentration de P700 dépend de la lumière et de la stabilité
        synthesis = self.synthesis_rate * light_intensity * self.stability
        degradation = (1 - self.stability) * self.p700_concentration * 0.1
        self.p700_concentration += (synthesis - degradation) * dt
        self.p700_concentration = np.clip(self.p700_concentration, 0, 1)

        # Le facteur de cohérence est lié à la stabilité
        self.coherence_factor = self.stability
