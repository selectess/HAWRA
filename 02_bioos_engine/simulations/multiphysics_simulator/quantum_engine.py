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
        
        # État quantique simplifié (amplitude de probabilité pour |0> et |1>)
        self.qubit_state = [1.0, 0.0] # Initialement dans l'état |0>

        # Paramètres du modèle
        self.p700_excitation_threshold = config.get('p700_threshold', 0.8)
        self.decoherence_rate = config.get('decoherence_rate', 0.01)

    def apply_gate(self, gate_name, qubits):
        """Applique une porte quantique à l'état du qubit."""
        print(f"QUANTUM_OP: Applying {gate_name} to {qubits}")
        if gate_name == 'H':
            # Hadamard simple: |0> -> (|0>+|1>)/sqrt(2), |1> -> (|0>-|1>)/sqrt(2)
            # Ici on simule l'effet sur les probabilités de mesure
            if self.qubit_state == [1.0, 0.0] or self.qubit_state == [0.0, 1.0]:
                self.qubit_state = [0.5, 0.5] # Superposition équilibrée
            else:
                self.qubit_state = [1.0, 0.0] # Retour à la base (simplifié)
        elif gate_name == 'X':
            self.qubit_state = [self.qubit_state[1], self.qubit_state[0]]
        elif gate_name == 'CCNOT':
            # Toffoli simplified for current single-qubit focus
            # In a real multi-qubit engine, this would affect 3 qubits
            # For HAWRA v1, we treat it as a conditional flip if the "system" allows it
            print("QUANTUM_OP: CCNOT (Toffoli) - multi-qubit gates are experimental in HAWRA v1")
            self.qubit_state = [self.qubit_state[1], self.qubit_state[0]]
        elif gate_name == 'DJ_ORACLE_CONSTANT':
            # Oracle constant: f(x) = 0 or f(x) = 1. No change to the query qubit.
            print("QUANTUM_OP: DJ Oracle (Constant)")
            pass 
        elif gate_name == 'DJ_ORACLE_BALANCED':
            # Oracle balanced: f(0)=0, f(1)=1 (CNOT)
            print("QUANTUM_OP: DJ Oracle (Balanced)")
            # Simule l'effet d'un CNOT sur la phase du qubit cible
            self.qubit_state = [self.qubit_state[1], self.qubit_state[0]]

    def measure(self, qubit):
        """Mesure le qubit et déclenche la bioluminescence."""
        print(f"MEASURE: Qubit {qubit}")
        # La mesure réduit l'état
        if random.random() < self.qubit_state[0]:
            self.p700_state = 0
            self.luc_green_output = 1.0
            self.qubit_state = [1.0, 0.0]
        else:
            self.p700_state = 1
            self.luc_red_output = 1.0
            self.qubit_state = [0.0, 1.0]

    def update_state(self, p700_concentration):
        """Met à jour l'état de P700 en fonction de la concentration et de la décohérence."""
        # Réinitialisation des sorties de lecture (impulsions brèves)
        self.luc_green_output = max(0, self.luc_green_output - 0.2)
        self.luc_red_output = max(0, self.luc_red_output - 0.2)

        # Si le P700 est bas, la cohérence est perdue
        if p700_concentration < self.p700_excitation_threshold:
            # Décohérence vers l'état fondamental
            if self.qubit_state[1] > 0:
                decay = self.decoherence_rate * (1.0 - p700_concentration)
                self.qubit_state[1] = max(0, self.qubit_state[1] - decay)
                self.qubit_state[0] = 1.0 - self.qubit_state[1]

    def get_state(self):
        """Retourne l'état actuel pour le logging."""
        return {
            "p700_state": self.p700_state,
            "luc_green_output": self.luc_green_output,
            "luc_red_output": self.luc_red_output,
            "qubit_state": self.qubit_state
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
