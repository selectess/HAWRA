import numpy as np
import qutip
from qutip import basis, sigmax, sigmay, sigmaz, sigmam, mesolve, fidelity, Qobj
from qutip.qip.operations import snot
import numpy as np
from qutip import expect, basis

class QuantumCore:
    """
    Quantum engine for simulating the dynamics of the P700 bio-qubit.
    Integrates control physics (gates) and decoherence (noise).
    """
    def __init__(self, config):
        """
        Initializes the quantum engine with a configuration.
        
        Args:
            config (dict): Configuration dictionary containing physical
                           parameters like t1, t2, etc.
        """
        self.config = config
        self.t1 = config.get('t1', 1e-6)  # T1 relaxation time in seconds
        self.t2 = config.get('t2', 1e-6)  # T2 dephasing time in seconds
        self.silica_protection_factor = config.get('silica_protection_factor', 1.0)
        self.options = {'nsteps': 50000}

        # Base states
        self.psi0 = basis(2, 0)
        self.psi1 = basis(2, 1)

    def _get_hamiltonian(self, control_pulse=None):
        """Constructs the Hamiltonian for the simulation."""
        if control_pulse:
            # Hamiltonian for a Rabi pulse: H = (Omega/2) * sigma_x
            # Omega (Rabi freq) = 2*pi*Amplitude. So H = pi*A*sigma_x
            amplitude = control_pulse.get('amplitude', 0)
            return np.pi * amplitude * sigmax()
        else:
            # Null Hamiltonian for pure decoherence simulations
            return 0 * sigmaz()

    def _get_collapse_operators(self):
        """Constructs the list of collapse operators to model noise."""
        # Relaxation rate (T1)
        gamma_relax = 1.0 / self.t1
        c_ops = [np.sqrt(gamma_relax) * sigmam()] if self.t1 > 0 else []

        # Dephasing rate (T2)
        # The total dephasing rate is 1/T2 = 1/(2*T1) + 1/T2*
        # We model T2* (pure dephasing) directly
        if self.t2 > 0:
            gamma_phi = (1.0 / self.t2) * self.silica_protection_factor
            c_ops.append(np.sqrt(gamma_phi) * sigmaz())
            
        return c_ops

    def get_gate(self, gate_name: str) -> Qobj:
        """Returns the matrix of the specified quantum gate."""
        if gate_name == 'H':
            return snot()
        elif gate_name == 'X':
            return sigmax()
        elif gate_name == 'Y':
            return sigmay()
        elif gate_name == 'Z':
            return sigmaz()
        elif gate_name == 'CNOT':
            # Approximation for a CNOT on a register with one qubit (placeholder)
            # We apply an action equivalent to X to illustrate the effect.
            return sigmax()
        else:
            raise ValueError(f"Unknown quantum gate: {gate_name}")

    def run_simulation(self, initial_state, duration, control_pulse=None, observables=None):
        """
        Runs the quantum dynamics simulation.

        Args:
            initial_state (Qobj): The initial quantum state.
            duration (float): The duration of the simulation in seconds.
            control_pulse (dict, optional): Describes the control pulse.
            observables (list, optional): List of operators for which to measure the expected value.

        Returns:
            qutip.Result: The result object from mesolve.
        """
        hamiltonian = self._get_hamiltonian(control_pulse)
        c_ops = self._get_collapse_operators()
        # Reduction of the internal resolution to speed up the simulation
        times = np.linspace(0, duration, 100)

        if observables is None:
            observables = []

        result = mesolve(hamiltonian, initial_state, times, c_ops=c_ops, e_ops=observables, options=self.options)
        return result

    def measure_qubit(self, qubit_state, basis='Z'):
        """
        Measures the qubit in a given basis.

        Args:
            qubit_state (Qobj): The current state of the qubit.
            basis (str): The measurement basis ('X', 'Y', or 'Z').

        Returns:
            (Qobj, int): The new state of the qubit after the measurement and the result (0 or 1).
        """
        if basis == 'Z':
            P0 = qutip.basis(2, 0) * qutip.basis(2, 0).dag()
            P1 = qutip.basis(2, 1) * qutip.basis(2, 1).dag()
        elif basis == 'X':
            P0 = (qutip.basis(2, 0) + qutip.basis(2, 1)).unit() * (qutip.basis(2, 0) + qutip.basis(2, 1)).unit().dag()
            P1 = (qutip.basis(2, 0) - qutip.basis(2, 1)).unit() * (qutip.basis(2, 0) - qutip.basis(2, 1)).unit().dag()
        elif basis == 'Y':
            P0 = (basis(2, 0) + 1j * basis(2, 1)).unit() * (basis(2, 0) + 1j * basis(2, 1)).unit().dag()
            P1 = (basis(2, 0) - 1j * basis(2, 1)).unit() * (basis(2, 0) - 1j * basis(2, 1)).unit().dag()
        else:
            raise ValueError("The measurement basis must be 'X', 'Y', or 'Z'.")

        prob0 = float(np.real(expect(P0, qubit_state)))
        prob1 = float(np.real(expect(P1, qubit_state)))

        if np.random.rand() < prob0:
            result = 0
            new_state = (P0 * qubit_state).unit()
        else:
            result = 1
            new_state = (P1 * qubit_state).unit()

        return new_state, result


class QuantumState:
    """
    Represents the state of a multi-qubit quantum system.

    This class manages the state vector of the quantum system and provides
    methods for applying gates, measuring qubits, and applying noise.
    """
    def __init__(self, qubit_count):
        """
        Initializes the quantum state.

        Args:
            qubit_count (int): The number of qubits in the system.
        """
        self.qubit_count = qubit_count
        self.state_vector = np.zeros(2**self.qubit_count, dtype=complex)
        self.state_vector[0] = 1.0

    def get_state(self):
        """
        Returns the current state of the quantum system.

        Returns:
            list: A list of dictionaries representing the complex amplitudes of
                the state vector.
        """
        return [{"real": c.real, "imag": c.imag} for c in self.state_vector]

    def _get_operator(self, gate_matrix, target_qubits):
        """
        Constructs the full operator for a given gate and target qubits.

        Args:
            gate_matrix (np.ndarray): The matrix of the gate to apply.
            target_qubits (list): A list of the target qubit indices.

        Returns:
            np.ndarray: The full operator for the given gate and target qubits.
        """
        num_target_qubits = len(target_qubits)
        op_list = [np.identity(2) for _ in range(self.qubit_count)]

        if num_target_qubits == 1:
            op_list[target_qubits[0]] = gate_matrix
        elif num_target_qubits == 2: # Specific for CNOT
            control, target = target_qubits[0], target_qubits[1]
            # This is a simplified construction and assumes CNOT. A more general method is needed for other 2-qubit gates.
            # We will build the operator that acts on the whole state space.
            I = np.identity(2)
            X = np.array([[0, 1], [1, 0]])
            P0 = np.array([[1, 0], [0, 0]]) # Projector |0><0|
            P1 = np.array([[0, 0], [0, 1]]) # Projector |1><1|

            # Operator for the control qubit part
            term1_ops = [I]*self.qubit_count
            term1_ops[control] = P0

            # Operator for the target qubit part
            term2_ops = [I]*self.qubit_count
            term2_ops[control] = P1
            term2_ops[target] = X

            # Construct the full operators via tensor product
            op1 = term1_ops[0]
            for i in range(1, self.qubit_count):
                op1 = np.kron(op1, term1_ops[i])

            op2 = term2_ops[0]
            for i in range(1, self.qubit_count):
                op2 = np.kron(op2, term2_ops[i])

            return op1 + op2

        # Construct the full operator via tensor product
        full_op = op_list[0]
        for i in range(1, self.qubit_count):
            full_op = np.kron(full_op, op_list[i])
        return full_op

    def apply_gate(self, gate_matrix, target_qubits):
        """
        Applies a quantum gate to the state vector.

        Args:
            gate_matrix (np.ndarray): The matrix of the gate to apply.
            target_qubits (list): A list of the target qubit indices.
        """
        num_target_qubits = len(target_qubits)
        
        # For CNOT, the gate_matrix passed is 4x4, but we need to build the full operator.
        # The logic is now inside _get_operator for 2-qubit gates.
        if num_target_qubits == 2:
            # We pass the target qubit indices, not the matrix itself, to the specialized constructor
            full_operator = self._get_operator(None, target_qubits)
        else:
            # For 1-qubit gates
            full_operator = self._get_operator(gate_matrix, target_qubits)

        self.state_vector = full_operator @ self.state_vector

    def measure(self, target_qubits):
        """
        Measures a qubit in the Z-basis.

        Args:
            target_qubits (list): A list of the target qubit indices.

        Returns:
            int: The measurement outcome (0 or 1).
        """
        # Simplified measurement: measures the first qubit in the list in the Z-basis.
        # A full implementation would handle multi-qubit measurements and different bases.
        qubit_to_measure = target_qubits[0]
        
        # Probability of measuring |0>
        prob0 = 0
        for i in range(len(self.state_vector)):
            # If the bit corresponding to qubit_to_measure is 0
            if (i >> (self.qubit_count - 1 - qubit_to_measure)) & 1 == 0:
                prob0 += np.abs(self.state_vector[i])**2

        # Randomly choose outcome based on probability
        if np.random.rand() < prob0:
            outcome = 0
            norm = np.sqrt(prob0)
            # Project state to |0>
            for i in range(len(self.state_vector)):
                if (i >> (self.qubit_count - 1 - qubit_to_measure)) & 1 != 0:
                    self.state_vector[i] = 0
        else:
            outcome = 1
            norm = np.sqrt(1 - prob0)
            # Project state to |1>
            for i in range(len(self.state_vector)):
                if (i >> (self.qubit_count - 1 - qubit_to_measure)) & 1 == 0:
                    self.state_vector[i] = 0

        # Normalize the state vector
        self.state_vector /= norm
        return outcome

    def apply_noise(self, coherence_factor):
        """Applies phase damping noise based on the coherence factor."""
        if coherence_factor >= 1.0:
            return # No noise if coherence is perfect

        # The strength of the noise is inversely related to the coherence factor
        noise_strength = (1.0 - coherence_factor) * 0.1 # Max noise strength of 0.1 radians

        for i in range(self.qubit_count):
            # Apply a random phase kick around the Z-axis to each qubit
            random_phase = np.exp(1j * np.random.uniform(-noise_strength, noise_strength))
            phase_op = np.array([[1, 0], [0, random_phase]])
            
            # Construct the full operator for this single qubit phase kick
            full_operator = self._get_operator(phase_op, [i])
            self.state_vector = full_operator @ self.state_vector
