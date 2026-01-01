# Quantum Engine for the Multiphysics Simulator

import numpy as np
import random

class QuantumEngine:
    """
    Simule la dynamique du cœur quantique P700 et l'architecture de lecture à double canal.
    Basé sur le design génétique de HAWRA_FINAL_VALIDATED.gb (v1).
    """
    def __init__(self, config):
        print("Initializing Quantum Engine with dual-channel readout model...")
        # État de P700: 0=fondamental, 1=excité (P700*)
        self.p700_state = 0
        # Sorties des canaux de lecture
        self.luc_green_output = 0.0 # Canal |0> (stable)
        self.luc_red_output = 0.0   # Canal |1> (instable)

        # Paramètres du modèle
        self.p700_excitation_threshold = config.get('p700_threshold', 0.8)
        self.decoherence_rate_slow = config.get('decoherence_rate_slow', 0.1) # Effondrement lent -> ATP
        self.decoherence_rate_fast = config.get('decoherence_rate_fast', 0.8) # Effondrement rapide -> ROS
        self.fast_collapse_probability = config.get('fast_collapse_probability', 0.2) # Probabilité d'un effondrement rapide

    def update_state(self, p700_concentration):
        """Met à jour l'état de P700 et des canaux de lecture."""
        # Réinitialisation des sorties à chaque pas
        self.luc_green_output = 0.0
        self.luc_red_output = 0.0

        if self.p700_state == 0:
            # Tente d'exciter P700 si la concentration est suffisante
            if p700_concentration > self.p700_excitation_threshold:
                excitation_prob = (p700_concentration - self.p700_excitation_threshold) / (1.0 - self.p700_excitation_threshold)
                if random.random() < excitation_prob:
                    self.p700_state = 1
                    # print("P700 EXCITÉ")
        else: # p700_state == 1
            # P700* est excité, il peut s'effondrer
            if random.random() < self.fast_collapse_probability:
                # Effondrement RAPIDE (simule production de ROS)
                if random.random() < self.decoherence_rate_fast:
                    self.p700_state = 0
                    self.luc_red_output = 1.0 # Le canal rouge s'active
                    # print("EFFONDREMENT RAPIDE -> CANAL ROUGE")
            else:
                # Effondrement LENT (simule production d'ATP)
                if random.random() < self.decoherence_rate_slow:
                    self.p700_state = 0
                    self.luc_green_output = 1.0 # Le canal vert s'active
                    # print("EFFONDREMENT LENT -> CANAL VERT")

    def get_state(self):
        """Retourne l'état actuel pour le logging."""
        return {
            "p700_state": self.p700_state,
            "luc_green_output": self.luc_green_output,
            "luc_red_output": self.luc_red_output
        }

    def update(self, time, dt, bio_state):
        p700_concentration = bio_state.get('p700_concentration', 0)

        # Probability of excitation is proportional to P700 concentration
        excitation_prob = np.clip(p700_concentration / self.threshold, 0, 1) * dt

        # Probability of decay (decoherence)
        decay_prob = self.decoherence_rate * dt

        # Get current probabilities
        p_ground, p_excited = self.qubit_state

        # State transitions based on probabilities
        if np.random.rand() < excitation_prob * p_ground:
            # Transition from ground to excited
            transition_amount = min(p_ground, 0.1) # Limit transition amount per step
            self.qubit_state = [p_ground - transition_amount, p_excited + transition_amount]
            print(f"EVENT: Qubit excitation with probability {excitation_prob:.2f}")
        elif np.random.rand() < decay_prob * p_excited:
            # Transition from excited to ground
            transition_amount = min(p_excited, 0.1) # Limit transition amount per step
            self.qubit_state = [p_ground + transition_amount, p_excited - transition_amount]
            print(f"EVENT: Qubit decay with probability {decay_prob:.2f}")

        # Normalize to ensure it remains a valid probability distribution
        norm = sum(self.qubit_state)
        if norm > 0:
            self.qubit_state = [p / norm for p in self.qubit_state]

        print(f"Updating quantum state at t={time}. Qubit state probabilities: P(G)={self.qubit_state[0]:.2f}, P(E)={self.qubit_state[1]:.2f}")
        
        return {
            'qubit_state': self.qubit_state
        }
