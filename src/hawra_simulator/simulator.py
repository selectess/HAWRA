import json
import numpy as np
from .engines.quantum_core import QuantumState, QuantumCore
from .engines.biological_system import BiologicalSystem
from .engines.environment import Environment

class Simulator:
    """
    The main simulator for the HAWRA project.

    This class is responsible for loading a BSIM script, running the simulation,
    and saving the results.
    """
    def __init__(self, bsim_script=None):
        """
        Initializes the simulator.

        Args:
            bsim_script (str, optional): The path to the BSIM script.
                Defaults to None.
        """
        self.bsim_script = bsim_script
        self.qubit_count = 0
        self.qubits = []
        self.quantum_state = QuantumState(self.qubit_count)
        self.results = []
        self.final_state = None
        self.time = 0.0
        self.dt = 1.0  # Default time step
        self.duration = 100.0  # Default simulation time

        self.environment = Environment({})
        self.biological_system = BiologicalSystem({})

        self.environment_history = []
        self.biological_system_history = []
        self.quantum_history = []
        self.quantum_core = None
        self._pulses = []

        self.gate_matrices = {
            'H': np.array([[1, 1], [1, -1]]) / np.sqrt(2),
            'X': np.array([[0, 1], [1, 0]]),
            'Y': np.array([[0, -1j], [1j, 0]]),
            'Z': np.array([[1, 0], [0, -1]]),
            'CNOT': np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
        }

        self.instructions = []
        if isinstance(bsim_script, dict):
            cfg = bsim_script
            env_cfg = cfg.get('environment', {})
            bio_cfg = cfg.get('biological_system', {})
            self.environment = Environment(env_cfg)
            self.biological_system = BiologicalSystem(bio_cfg)
        elif self.bsim_script:
            self.load_bsim_script(self.bsim_script)

    def load_bsim_script(self, script_path):
        """
        Loads a BSIM script from a JSON file.

        Args:
            script_path (str): The path to the BSIM script.
        """
        with open(script_path, 'r') as f:
            script = json.load(f)
        
        meta = script.get('metadata', {})
        self.qubit_count = int(meta.get('qubit_count', 1))
        self.qubits = meta.get('qubits', [])
        self.instructions = script['instructions']
        # Normalize alternative instruction schemas
        norm = []
        for instr in self.instructions:
            if 'command' in instr:
                norm.append(instr)
                continue
            op = instr.get('op_code') or instr.get('opcode')
            p = instr.get('params', {})
            if not op:
                norm.append({'command': None, 'params': p})
                continue
            if op == 'APPLY_GATE':
                gate = p.get('gate')
                targets = p.get('target') or p.get('targets') or []
                norm.append({'command': 'QUANTUM_OP', 'params': {'gate': gate, 'qubits': targets}})
            elif op == 'MEASURE':
                targets = p.get('target') or p.get('targets') or []
                cbit = p.get('cbit') or p.get('classical_bit')
                qubit = targets[0] if isinstance(targets, list) and targets else targets
                norm.append({'command': 'MEASURE', 'params': {'qubit': qubit, 'classical_bit': cbit[0] if isinstance(cbit, list) and cbit else cbit}})
            else:
                norm.append({'command': op, 'params': p})
        self.instructions = norm
        try:
            self.quantum_state = QuantumState(self.qubit_count)
        except Exception:
            self.quantum_state = None

        # Initialize engines from config
        sim_config = script.get('simulation_config', {})
        self.duration = sim_config.get('duration', 100.0)
        self.dt = sim_config.get('dt', 1.0)
        self.environment_script = script.get('environment_script', [])
        
        env_config = script.get('environment_config', {})
        bio_config = script.get('biological_config', {})
        quantum_config = script.get('quantum_config', {})
        
        self.environment = Environment(env_config)
        self.biological_system = BiologicalSystem(bio_config)
        self.quantum_core = QuantumCore(quantum_config or {})

    def run_script(self, instructions):
        """
        Runs a BSIM instruction list directly (INITIALIZE + steps/stimuli),
        preparing environment and biological configs from the INITIALIZE block.
        """
        if not instructions:
            return {}
        init = instructions[0]
        if init.get('command') != 'INITIALIZE':
            # Fallback: assume defaults
            self.dt = 1.0
            self.environment = Environment({})
            self.biological_system = BiologicalSystem({'genes': [], 'grn': {}})
        else:
            cfg = init.get('config', {})
            self.dt = float(cfg.get('dt', 1.0))
            # Prepare environment
            env_cfg = cfg.get('env', {})
            self.environment = Environment(env_cfg)
            # Prepare biological config, adapting structures
            bio_cfg = cfg.get('bio', {})
            quantum_cfg = cfg.get('quantum', {})
            genes_dict = bio_cfg.get('genes', {})
            genes_list = []
            for gid, gparams in genes_dict.items():
                gentry = {'id': gid}
                gentry.update(gparams)
                gentry.setdefault('initial_expression', 0.0)
                genes_list.append(gentry)
            grn_cfg = bio_cfg.get('grn', {})
            grn_adapted = {}
            for tgt, info in grn_cfg.items():
                regs = []
                for reg in info.get('regulators', []):
                    regs.append({
                        'id': reg.get('gene'),
                        'type': reg.get('type'),
                        'weight': reg.get('weight', 1.0),
                        'hill_coefficient': reg.get('hill_coefficient', 1.0),
                        'half_max_concentration': reg.get('half_max_concentration', 0.5)
                    })
                grn_adapted[tgt] = regs
            bio_init = {
                'genes': genes_list,
                'grn': grn_adapted,
                'synthesis_rate': bio_cfg.get('p700_synthesis_rate', 0.05),
                'initial_p700': 1.0
            }
            self.biological_system = BiologicalSystem(bio_init)
            self.quantum_core = QuantumCore(quantum_cfg or {})
        # Assign instructions and run
        self.instructions = instructions
        return self.run()


    def _execute_environment_script(self, time, dt):
        """
        Executes the environment script events that occur in the current time step.

        Args:
            time (float): The current simulation time.
            dt (float): The time step.
        """
        if not self.environment_script:
            return

        end_time = time + dt
        for event in self.environment_script:
            event_time = event['time']
            if time < event_time <= end_time:
                for action, params in event['actions'].items():
                    if hasattr(self.environment, action):
                        method_to_call = getattr(self.environment, action)
                        method_to_call(params)
                    else:
                        print(f"Avertissement : L'action '{action}' n'est pas définie dans l'environnement.")

    def run(self, duration=None, dt=None):
        """
        Runs the simulation.

        Returns:
            dict: A dictionary containing the simulation results.
        """
        if self.instructions:
            print("Début de l'exécution du script BSIM...")
            self.time = 0.0
            for instruction in self.instructions:
                command = instruction.get('command')
                params = instruction.get('params', {})
                
                if command == 'RUN_UNTIL':
                    target_time = params.get('time')
                    print(f"[RUN] Exécution jusqu'au temps {target_time}...")
                    
                    for t in np.arange(self.time, target_time, self.dt):
                        self.environment.update_from_schedule(t)
                        env_state = self.environment.get_state()
                        self.environment_history.append({'time': t, **env_state})

                        self.biological_system.update(self.dt, env_state['light_intensity'], env_state.get('light_wavelength'))
                        bio_state = self.biological_system.get_state()
                        self.biological_system_history.append({
                            'time': t,
                            'light_intensity': env_state.get('light_intensity', 0.0),
                            'coherence_factor': bio_state['coherence_factor'],
                            'p700_concentration': bio_state['p700_concentration']
                        })
                    
                    self.time = target_time

                elif command == 'STIMULUS_APPLY':
                    print(f"[STIMULUS] Application d'un stimulus: {params}")
                    stim = params.get('stimulus')
                    if stim == 'light_pulse':
                        args = params.get('arguments', {})
                        intensity = float(args.get('intensity', 0.0))
                        duration = float(args.get('duration', 0.0))
                        start_time = self.time
                        if duration <= 0.0 or intensity < 0.0:
                            raise ValueError("Invalid light_pulse parameters: duration must be > 0 and intensity >= 0")
                        end_time = start_time + duration
                        self._pulses.append({'start': start_time, 'end': end_time, 'intensity': intensity})
                        times = sorted({0.0, *[p['start'] for p in self._pulses], *[p['end'] for p in self._pulses]})
                        schedule = []
                        current = None
                        for t in times:
                            active = [p for p in self._pulses if p['start'] <= t < p['end']]
                            if active:
                                idx = max(range(len(active)), key=lambda i: active[i]['start'])
                                val = active[idx]['intensity']
                                wav = active[idx].get('wavelength')
                            else:
                                val = 0.0
                                wav = None
                            if current is None or val != current:
                                entry = {'time': t, 'intensity': val}
                                if wav is not None:
                                    entry['wavelength'] = float(wav)
                                schedule.append(entry)
                                current = val
                        self.environment.light_schedule = schedule
                        if not hasattr(self.environment, 'light_schedule') or self.environment.light_schedule is None:
                            self.environment.light_schedule = []
                        if hasattr(self.environment, 'compact_schedule'):
                            self.environment.compact_schedule()
                    elif stim == 'adaptive_light':
                        args = params.get('arguments', {})
                        intensity = float(args.get('intensity', 0.0))
                        duration = float(args.get('duration', 0.0))
                        wavelength = float(args.get('wavelength', 660.0))
                        start_time = self.time
                        if duration <= 0.0 or intensity < 0.0:
                            raise ValueError("Invalid adaptive_light parameters: duration must be > 0 and intensity >= 0")
                        end_time = start_time + duration
                        self._pulses.append({'start': start_time, 'end': end_time, 'intensity': intensity, 'wavelength': wavelength})
                        times = sorted({0.0, *[p['start'] for p in self._pulses], *[p['end'] for p in self._pulses]})
                        schedule = []
                        current = None
                        for t in times:
                            active = [p for p in self._pulses if p['start'] <= t < p['end']]
                            if active:
                                idx = max(range(len(active)), key=lambda i: active[i]['start'])
                                val = active[idx]['intensity']
                                wav = active[idx].get('wavelength')
                            else:
                                val = 0.0
                                wav = None
                            if current is None or val != current:
                                entry = {'time': t, 'intensity': val}
                                if wav is not None:
                                    entry['wavelength'] = float(wav)
                                schedule.append(entry)
                                current = val
                        self.environment.light_schedule = schedule
                        if not hasattr(self.environment, 'light_schedule') or self.environment.light_schedule is None:
                            self.environment.light_schedule = []
                        if hasattr(self.environment, 'compact_schedule'):
                            self.environment.compact_schedule()
                    elif stim == 'light':
                        args = params.get('arguments', {})
                        intensity = float(args.get('intensity', 0.0))
                        self.environment.set_light_intensity(intensity)
                        if 'wavelength' in args:
                            self.environment.set_light_wavelength(float(args.get('wavelength')))
                    elif stim == 'heat':
                        args = params.get('arguments', {})
                        temperature = float(args.get('temperature', self.environment.temperature))
                        self.environment.set_temperature(temperature)
                    else:
                        self.biological_system.apply_stimulus(params)

                elif command == 'QUANTUM_OP':
                    p = instruction.get('params', {})
                    gate = p.get('gate')
                    qubits = p.get('qubits', [])
                    if self.quantum_state is None:
                        continue
                    if isinstance(qubits, list) and qubits and isinstance(qubits[0], int):
                        targets = qubits
                    elif isinstance(qubits, list) and qubits:
                        name_to_index = {}
                        for idx, q in enumerate(self.qubits):
                            n = q.get('name') if isinstance(q, dict) else None
                            name_to_index[n] = idx
                        targets = [name_to_index.get(q, 0) for q in qubits]
                    else:
                        targets = [0] if self.qubit_count >= 1 else []
                    gm = self.gate_matrices.get(gate)
                    if gm is not None and targets:
                        self.quantum_state.apply_gate(gm, targets)
                        self.quantum_history.append({'time': self.time, 'op': gate, 'targets': targets})

                elif command == 'MEASURE':
                    p = instruction.get('params', {})
                    target = p.get('qubit')
                    if self.quantum_state is None:
                        continue
                    if isinstance(target, int):
                        targets = [target]
                    else:
                        name_to_index = {}
                        for idx, q in enumerate(self.qubits):
                            n = q.get('name') if isinstance(q, dict) else None
                            name_to_index[n] = idx
                        targets = [name_to_index.get(target, 0)]
                    outcome = self.quantum_state.measure(targets)
                    self.quantum_history.append({'time': self.time, 'measure': targets[0], 'outcome': int(outcome)})

                elif command == 'run':
                    circ = instruction.get('circuit')
                    args = instruction.get('arguments', {})
                    amp = float(args.get('amplitude', 0.0)) if isinstance(args.get('amplitude', 0.0), (int, float)) else 0.0
                    dur = float(args.get('duration', 0.01)) if isinstance(args.get('duration', 0.01), (int, float)) else 0.01
                    if self.quantum_core:
                        from qutip import sigmaz
                        result = self.quantum_core.run_simulation(self.quantum_core.psi0, dur, control_pulse={'amplitude': amp}, observables=[sigmaz()])
                        z_raw = result.expect[0] if hasattr(result, 'expect') and result.expect else []
                        z_arr = np.asarray(z_raw).reshape(-1)
                        z_series = [float(x) for x in z_arr.tolist()]
                        self.quantum_history.append({'time': self.time, 'run': circ, 'args': args, 'z_exp': z_series})
                    else:
                        self.quantum_history.append({'time': self.time, 'run': circ, 'args': args})

                elif command == 'INITIALIZE':
                    print("[INIT] Configuration du simulateur...")
                    pass

                else:
                    print(f"[WARN] Instruction non reconnue ou non implémentée: {command}")

            final_state = self.biological_system.get_state()
            results = {
                'results': self.biological_system_history,
                'final_state': final_state,
                'environment_history': self.environment_history,
                'quantum_history': self.quantum_history
            }
            return results

        use_duration = float(duration) if duration is not None else float(self.duration)
        use_dt = float(dt) if dt is not None else float(self.dt)
        env_hist = []
        bio_hist = []
        t = 0.0
        steps = int(use_duration / use_dt)
        for _ in range(steps):
            if hasattr(self.environment, 'update_from_schedule'):
                self.environment.update_from_schedule(t)
            env_state = self.environment.get_state()
            env_hist.append({'time': t, **env_state})
            self.biological_system.update(use_dt, env_state['light_intensity'], env_state.get('light_wavelength'))
            bio_state = self.biological_system.get_state()
            bio_hist.append({
                'time': t,
                'light_intensity': env_state.get('light_intensity', 0.0),
                'coherence_factor': bio_state['coherence_factor'],
                'p700_concentration': bio_state['p700_concentration']
            })
            t += use_dt
        final_state = self.biological_system.get_state()
        return {
            'results': bio_hist,
            'final_state': final_state,
            'environment_history': env_hist,
            'quantum_history': self.quantum_history
        }

    def save_results(self, output_path, simulation_results):
        """
        Saves the simulation results to a JSON file.

        Args:
            output_path (str): The path to the output file.
            simulation_results (dict): The simulation results.
        """
        def serialize_complex(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, complex):
                return {'real': obj.real, 'imag': obj.imag}
            return obj

        with open(output_path, 'w') as f:
            json.dump(simulation_results, f, indent=4, default=serialize_complex)
        print(f"Résultats de la simulation sauvegardés dans {output_path}")
