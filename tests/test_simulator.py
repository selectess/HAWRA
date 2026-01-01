import unittest
import numpy as np
from qutip import basis, sigmax, fidelity
import os
import sys
import importlib
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src', 'simulator')))

QuantumCore = importlib.import_module('hawra_simulator.engines.quantum_core').QuantumCore
BiologicalSystem = importlib.import_module('hawra_simulator.engines.biological_system').BiologicalSystem
Simulator = importlib.import_module('hawra_simulator.simulator').Simulator

class TestQuantumCore(unittest.TestCase):

    def setUp(self):
        """Configuration initiale pour les tests du QuantumCore."""
        self.config = {
            't1': 1.5e-9,  # 1.5 ns
            't2': 0.5e-9,  # 0.5 ns
            'silica_protection_factor': 1/3 # Facteur d'amélioration de T2 (T2_new = T2_old * 3)
        }
        self.quantum_engine = QuantumCore(self.config)

    def test_decoherence_simulation(self):
        """Teste la simulation de décohérence pure (pas de contrôle)."""
        # État initial en superposition
        initial_state = (basis(2, 0) + basis(2, 1)).unit()
        duration = 2e-9 # 2 ns

        # Exécute la simulation sans impulsion de contrôle
        result = self.quantum_engine.run_simulation(initial_state, duration)

        # L'état final devrait avoir perdu sa cohérence et tendu vers un état mixte
        final_state = result.states[-1]
        
        # La fidélité avec l'état initial doit être faible
        f = fidelity(initial_state, final_state)
        self.assertLess(f, 0.8, "La cohérence devrait chuter de manière significative.")

    def test_pi_pulse_simulation(self):
        """Teste l'application d'une impulsion Pi (porte NOT)."""
        # État initial |0>
        initial_state = basis(2, 0)
        
        # Une impulsion Pi doit inverser le qubit. Durée = pi / (pi * Amplitude)
        # Pour une porte NOT parfaite, Amplitude = 0.5/durée
        duration = 0.1e-9 # 100 ps, assez court pour minimiser la décohérence
        amplitude = 0.5 / duration
        control_pulse = {'amplitude': amplitude}

        # État final idéal |1>
        target_state = basis(2, 1)

        # Exécute la simulation
        result = self.quantum_engine.run_simulation(initial_state, duration, control_pulse)
        final_state = result.states[-1]

        # La fidélité avec l'état |1> doit être élevée
        f = fidelity(target_state, final_state)
        self.assertGreater(f, 0.95, "La fidélité de la porte NOT doit être élevée.")

class TestBiologicalSystem(unittest.TestCase):

    def setUp(self):
        """Configuration initiale pour les tests du BiologicalSystem."""
        self.config = {'optimal_temp': 35.0}
        self.bio_system = BiologicalSystem(self.config)

    def test_temperature_stability(self):
        """Vérifie que la stabilité est maximale à la température optimale."""
        # À la température optimale
        state_optimal = self.bio_system.update_state(35.0, dt=1.0)
        self.assertAlmostEqual(state_optimal['stability'], 1.0, places=5, msg="La stabilité doit être maximale à la température optimale.")
        self.assertAlmostEqual(state_optimal['coherence_factor'], 1.0, places=5, msg="Le facteur de cohérence doit être maximal à la température optimale.")

        # En dehors de la température optimale
        state_off_optimal = self.bio_system.update_state(40.0, dt=1.0)
        self.assertLess(state_off_optimal['stability'], 1.0, "La stabilité doit diminuer en dehors de la température optimale.")
        self.assertLess(state_off_optimal['coherence_factor'], 1.0, "Le facteur de cohérence doit diminuer en dehors de la température optimale.")

    def test_chemical_stimulus_up_down(self):
        cfg = {
            'genes': [
                {'id': 'x', 'basal_rate': 0.1, 'degradation_rate': 0.05, 'initial_expression': 0.0}
            ],
            'grn': {}
        }
        bio = BiologicalSystem(cfg)
        base = bio.gene_params['x']['basal_rate']
        bio.apply_stimulus({'stimulus': 'chemical', 'target': 'x', 'arguments': {'effect': 'upregulate', 'magnitude': 0.2}})
        self.assertGreater(bio.gene_params['x']['basal_rate'], base)
        base2 = bio.gene_params['x']['basal_rate']
        bio.apply_stimulus({'stimulus': 'chemical', 'target': 'x', 'arguments': {'effect': 'downregulate', 'magnitude': 0.5}})
        self.assertLess(bio.gene_params['x']['basal_rate'], base2)

class TestSimulator(unittest.TestCase):

    def test_full_simulation_run(self):
        """Test d'une exécution complète de la simulation de bout en bout."""
        config = {
            'quantum_core': {'t1': 50e-6, 't2': 20e-6},
            'biological_system': {'p700_synthesis_rate': 0.1},
            'environment': {'temperature': 35.0}
        }
        simulator = Simulator(config)
        
        duration = 1e-7
        dt = 1e-8
        output = simulator.run(duration=duration, dt=dt)
        log = output['results']

        # Vérifier que le log a le bon nombre d'entrées
        expected_steps = int(duration / dt)
        self.assertEqual(len(log), expected_steps)

        env_hist = output['environment_history']
        for e in env_hist:
            self.assertEqual(e['temperature'], 35.0)
        for entry in log:
            self.assertAlmostEqual(entry['coherence_factor'], 1.0, places=5)
        last_entry = log[-1]
        self.assertTrue(0 <= last_entry['p700_concentration'] <= 1)

    def test_spectral_pathway_response(self):
        instructions = [
            {
                'command': 'INITIALIZE',
                'config': {
                    'dt': 1.0,
                    'env': {},
                    'bio': {
                        'genes': {
                            'gA': {'basal_rate': 0.02, 'initial_expression': 0.05, 'light_sensitivity': 0.5, 'peak_wavelength': 650, 'peak_sigma': 25},
                            'gB': {'basal_rate': 0.02, 'initial_expression': 0.05, 'light_sensitivity': 0.5, 'peak_wavelength': 450, 'peak_sigma': 25}
                        },
                        'grn': {}
                    }
                }
            },
            {
                'command': 'STIMULUS_APPLY',
                'params': {'stimulus': 'adaptive_light', 'arguments': {'intensity': 1.0, 'duration': 10, 'wavelength': 650}}
            },
            {'command': 'RUN_UNTIL', 'params': {'time': 12}}
        ]
        sim = Simulator({})
        out = sim.run_script(instructions)
        final_genes = out['final_state']['genes']
        self.assertGreater(final_genes['gA']['expression_level'], final_genes['gB']['expression_level'])

    def test_bsim_compatibility_op_code(self):
        import os
        import sys
        sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
        from hawra_simulator.simulator import Simulator
        bsim_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'arbol', 'phytoqmmml_demo.bsim.json'))
        sim = Simulator(bsim_path)
        out = sim.run()
        qh = out.get('quantum_history', [])
        self.assertTrue(len(qh) > 0)

    def test_pulse_schedule_metrics(self):
        """Vérifie qu'un schedule de pulses produit des métriques raisonnables."""
        # Crée un schedule de pulses courts
        pulses = [(0, 0.0)]
        for t in range(10, 100, 10):
            pulses.append((t, 60.0))
            pulses.append((t+2, 0.0))
        config = {
            'environment': {'temperature': 35.0, 'light_schedule': pulses},
            'biological_system': {
                'initial_p700': 1.0,
                'genes': [
                    {'id': 'gene_light', 'basal_rate': 0.01, 'degradation_rate': 0.05, 'light_sensitivity': 0.5, 'initial_expression': 0.0},
                    {'id': 'gene_a', 'basal_rate': 0.01, 'degradation_rate': 0.05, 'light_sensitivity': 0.0, 'initial_expression': 0.0},
                    {'id': 'gene_b', 'basal_rate': 0.02, 'degradation_rate': 0.06, 'light_sensitivity': 0.0, 'initial_expression': 0.02}
                ],
                'grn': {
                    'gene_a': {'gene_light': 1.0},
                    'gene_b': {'gene_a': -1.0}
                },
                'synthesis_rate': 0.05
            }
        }
        sim = Simulator(config)
        output = sim.run(duration=100, dt=1)
        log = output['results']
        self.assertEqual(len(log), 100)
        # Moyenne de lumière ~ 10 à 12 pour pulses 2/10
        light_avg = sum(x.get('light_intensity', 0.0) for x in log) / len(log)
        self.assertGreater(light_avg, 8.0)
        self.assertLess(light_avg, 14.0)
        p700_avg = sum(x['p700_concentration'] for x in log) / len(log)
        self.assertGreater(p700_avg, 0.7)

    def test_stimulus_apply_updates_schedule(self):
        instructions = [
            {
                'instruction_id': 0,
                'command': 'INITIALIZE',
                'config': {
                    'dt': 1,
                    'env': {'light_schedule': [(0, 0.0)]},
                    'bio': {'genes': {}, 'grn': {}}
                }
            },
            {
                'instruction_id': 1,
                'command': 'STIMULUS_APPLY',
                'params': {
                    'stimulus': 'light_pulse',
                    'target': 'GENE_A',
                    'arguments': {'intensity': 10.0, 'duration': 5.0}
                }
            },
            {
                'instruction_id': 2,
                'command': 'RUN_UNTIL',
                'params': {'time': 6.0}
            }
        ]
        sim = Simulator({})
        output = sim.run_script(instructions)
        env_hist = output['environment_history']
        ints = [e['light_intensity'] for e in env_hist]
        self.assertTrue(all(v == 10.0 for v in ints[:5]))
        self.assertEqual(ints[5], 0.0)

    def test_metrics_parity_instructions_vs_schedule(self):
        instructions = [
            {
                'instruction_id': 0,
                'command': 'INITIALIZE',
                'config': {
                    'dt': 1,
                    'env': {'light_schedule': [(0, 0.0)]},
                    'bio': {'genes': {}, 'grn': {}, 'p700_synthesis_rate': 0.05}
                }
            },
            {
                'instruction_id': 1,
                'command': 'STIMULUS_APPLY',
                'params': {
                    'stimulus': 'light_pulse',
                    'target': 'GENE_A',
                    'arguments': {'intensity': 40.0, 'duration': 3.0}
                }
            },
            {
                'instruction_id': 2,
                'command': 'RUN_UNTIL',
                'params': {'time': 30.0}
            }
        ]
        sim_inst = Simulator({})
        out_inst = sim_inst.run_script(instructions)
        rec_inst = out_inst['results']
        sched = [(0, 0.0), (0, 40.0), (3, 0.0)]
        sim_sched = Simulator({'environment': {'temperature': 35.0, 'light_schedule': sched}, 'biological_system': {'initial_p700': 1.0}})
        out_sched = sim_sched.run(duration=30, dt=1)
        rec_sched = out_sched['results']
        def avg(records, key):
            return sum(r.get(key, 0.0) for r in records) / len(records)
        tol = 0.03
        self.assertLess(abs(avg(rec_inst, 'light_intensity') - avg(rec_sched, 'light_intensity'))/max(1e-9, avg(rec_sched, 'light_intensity')), tol)
        self.assertLess(abs(avg(rec_inst, 'p700_concentration') - avg(rec_sched, 'p700_concentration'))/max(1e-9, avg(rec_sched, 'p700_concentration')), 0.1)

if __name__ == '__main__':
    unittest.main()
