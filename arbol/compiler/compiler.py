import sys
import os
import json
from .lexer import Lexer
from .parser import Parser
from .arbol_ast import *
from .error import ErrorReporter, CompilationError

class Compiler:
    def __init__(self, error_reporter=None):
        self.assembly = {
            "metadata": {
                "source_arbol": "compiled.arbol",
                "version": "0.1"
            },
            "instructions": []
        }
        self.instruction_id = 0
        self.symbol_table = {
            'genes': {},
            'stimuli': {},
            'circuits': {},
            'qubits': {},
            'gates': {},
            'classical_bits': {}
        }
        self.error_reporter = error_reporter if error_reporter else ErrorReporter()
        self.initialize_assembly()

    def compile(self, ast):
        """Parcourt l'AST et génère l'assembly quantique."""
        try:
            self.visit(ast)
        except CompilationError as e:
            self.error_reporter.report(e.message, e.line, e.column)
        return self.assembly

    def add_instruction(self, command, params):
        instruction = {
            "instruction_id": self.instruction_id,
            "command": command,
            "params": params
        }
        self.assembly['instructions'].append(instruction)
        self.instruction_id += 1

    def visit(self, node):
        method_name = 'visit_' + type(node).__name__
        visitor = getattr(self, method_name, self.generic_visit)
        return visitor(node)

    def generic_visit(self, node):
        if isinstance(node, list):
            for item in node:
                self.visit(item)
        elif isinstance(node, Node):
            for field, value in node.__dict__.items():
                if isinstance(value, (Node, list)):
                    self.visit(value)

    def initialize_assembly(self):
        # This could be expanded to include more default configurations
        init_instruction = {
            "instruction_id": self.instruction_id,
            "command": "INITIALIZE",
            "config": {
                "max_time": 500, # Default max_time
                "dt": 5,
                "env": {"light_schedule": []},
                "quantum": {"p700_threshold": 0.8, "decoherence_rate": 0.05}
            }
        }
        config = init_instruction['config']
        if 'bio' not in config:
            config['bio'] = {
                'p700_synthesis_rate': 1.0,
                'p700_degradation_rate': 0.1,
                'optimal_temp': 35.0,
                'temp_sensitivity': 0.1,
                'genes': {},
                'grn': {}
            }

        self.assembly['instructions'].append(init_instruction)
        self.instruction_id += 1

    def visit_GenesBlock(self, node):
        for gene_def in node.genes:
            self.visit(gene_def)

    def visit_GeneDefinition(self, node):
        gene_name = node.gene_name.name
        if gene_name in self.symbol_table['genes']:
            self.error_reporter.report(f"Gene '{gene_name}' already defined.", 0, 0) # Line/col info missing
        self.symbol_table['genes'][gene_name] = {'defined': True}
        
        # Also add to the bio config
        if 'genes' not in self.assembly['instructions'][0]['config']['bio']:
            self.assembly['instructions'][0]['config']['bio']['genes'] = {}
        
        gene_config = {
            "basal_rate": 0.0,
            "degradation_rate": 0.1,
            "light_sensitivity": 0.0
        }

        for prop in node.properties:
            gene_config[prop.name.name] = prop.value

        self.assembly['instructions'][0]['config']['bio']['genes'][gene_name] = gene_config

    def visit_GRNDefinition(self, node):
        if 'grn' not in self.assembly['instructions'][0]['config']['bio']:
            self.assembly['instructions'][0]['config']['bio']['grn'] = {}
        grn = self.assembly['instructions'][0]['config']['bio']['grn']

        for interaction in node.interactions:
            target_gene = interaction.target.name
            if target_gene not in grn:
                grn[target_gene] = {"regulators": []}

            regulator_gene = interaction.regulator.name
            interaction_type = interaction.interaction_type
            # Normalize parameter keys to be case-insensitive
            raw_params = {p.name.name: p.value for p in interaction.params}
            params_norm = {str(k).lower(): v for k, v in raw_params.items()}

            grn[target_gene]["regulators"].append({
                "gene": regulator_gene,
                "type": "activator" if interaction_type == 'ACTIVATES' else "repressor",
                "weight": float(params_norm.get('weight', 1.0)),
                "hill_coefficient": float(params_norm.get('hill_coefficient', 2.0)),
                "half_max_concentration": float(params_norm.get('half_max_concentration', 0.5))
            })

    def visit_RunBlock(self, node):
        self.generic_visit(node.instructions)

    def visit_StepInstruction(self, node):
        duration = self.parse_duration(node.duration.value)
        if duration is None:
            self.error_reporter.report(f"Invalid duration format '{node.duration.value}'.", node.line, node.column)
            duration = 0.0

        last_time = self.get_last_time()
        run_until_time = last_time + duration

        run_instruction = {
            "instruction_id": self.instruction_id,
            "command": "RUN_UNTIL",
            "params": {"time": run_until_time}
        }
        self.assembly['instructions'].append(run_instruction)
        self.instruction_id += 1

    def visit_CircuitDefinition(self, node):
        if node.name.name in self.symbol_table['circuits']:
            self.error_reporter.report(f"Circuit '{node.name.name}' already defined.", node.line, node.column)
        self.symbol_table['circuits'][node.name.name] = node

    def visit_StimulusDefinition(self, node):
        stimulus_name = node.name.name
        if stimulus_name in self.symbol_table['stimuli']:
            self.error_reporter.report(f"Stimulus '{stimulus_name}' already defined.", node.line, node.column)
            return

        parameters = {param.name: param.type for param in node.parameters}
        stim_type = parameters.get('type')
        self.symbol_table['stimuli'][stimulus_name] = {
            'type': stim_type if stim_type else 'stimulus',
            'parameters': parameters,
            'line': node.line,
            'column': node.column
        }

    def visit_StimulusApplication(self, node):
        stimulus_name = node.stimulus_name.name
        target = node.target.name if node.target else None
        arguments = {arg.name.name: float(arg.value) for arg in node.arguments}

        self.add_instruction('STIMULUS_APPLY', {
            "stimulus": stimulus_name,
            "target": target,
            "arguments": arguments
        })


    def visit_RunCircuit(self, node):
        circuit_name = node.circuit_name.name
        if circuit_name not in self.symbol_table['circuits']:
            self.error_reporter.report(f"Circuit '{circuit_name}' is not defined.", node.line, node.column)
            return

        circuit_def = self.symbol_table['circuits'][circuit_name]
        defined_params = {param.name: param.type for param in circuit_def.parameters}
        provided_args = {arg.name.name: arg.value for arg in node.arguments}

        for arg_name in provided_args:
            if arg_name not in defined_params:
                self.error_reporter.report(f"Unknown argument '{arg_name}' in run statement for circuit '{circuit_name}'.", node.line, node.column)

        for param_name in defined_params:
            if param_name not in provided_args:
                self.error_reporter.report(f"Missing argument '{param_name}' in run statement for circuit '{circuit_name}'.", node.line, node.column)

        run_instruction = {
            "instruction_id": self.instruction_id,
            "command": "run",
            "circuit": circuit_name,
            "arguments": provided_args
        }
        self.assembly['instructions'].append(run_instruction)
        self.instruction_id += 1

        for stmt in circuit_def.body:
            self.visit(stmt)

    def get_last_time(self):
        times = [instr.get('params', {}).get('time')
                 for instr in self.assembly['instructions']
                 if instr.get('command') == 'RUN_UNTIL']
        if not times:
            return 0
        return max(times)

    def parse_duration(self, duration_str):
        if isinstance(duration_str, (int, float)):
            return float(duration_str)
        s = str(duration_str).replace('"', '').strip()
        try:
            return float(s)
        except ValueError:
            pass
        if s.endswith('ms'):
            return float(s[:-2]) / 1000.0
        if s.endswith('s'):
            return float(s[:-1])
        if s.endswith('min'):
            return float(s[:-3]) * 60.0
        if s.endswith('h'):
            return float(s[:-1]) * 3600.0
        return None

    def visit_LogicalQubitDefinition(self, node):
        qubit_name = node.name.name
        if qubit_name in self.symbol_table['qubits']:
            self.error_reporter.report(f"Qubit '{qubit_name}' already defined.", node.line, node.column)
        self.symbol_table['qubits'][qubit_name] = {
            'is_clause': node.is_clause,
            'line': node.line,
            'column': node.column
        }

    def visit_GateDefinition(self, node):
        gate_key = (node.name.name, tuple(ref.name for ref in node.qubit_refs))
        if gate_key in self.symbol_table['gates']:
            self.error_reporter.report(f"Gate '{node.name.name}' on {node.qubit_refs} already defined.", node.line, node.column)
        self.symbol_table['gates'][gate_key] = node.body

    def visit_QuantumOperation(self, node):
        gate_key = (node.name.name, tuple(ref.name for ref in node.qubit_refs))
        if gate_key not in self.symbol_table['gates']:
            if node.name.name.upper() in ['H', 'X', 'Y', 'Z', 'CNOT']:
                op_instruction = {
                    "instruction_id": self.instruction_id,
                    "command": "QUANTUM_OP",
                    "params": {
                        "gate": node.name.name.upper(),
                        "qubits": [ref.name for ref in node.qubit_refs]
                    }
                }
                self.assembly['instructions'].append(op_instruction)
                self.instruction_id += 1
                return
            self.error_reporter.report(f"Gate '{node.name.name}' on {node.qubit_refs} is not defined.", node.line, node.column)
            return
        
        implementation = self.symbol_table['gates'][gate_key]
        self.visit(implementation)

    def visit_GateApplication(self, node):
        gate_name = node.gate_name.name
        targets = [t.name for t in node.qubits]
        gate_key = (gate_name, tuple(targets))

        gate_name_upper = gate_name.upper()
        if gate_name_upper in ['H', 'X', 'Y', 'Z', 'CNOT', 'CX']:
            gate_to_emit = 'CNOT' if gate_name_upper == 'CX' else gate_name_upper
            op_instruction = {
                "instruction_id": self.instruction_id,
                "command": "QUANTUM_OP",
                "params": {
                    "gate": gate_to_emit,
                    "qubits": targets
                }
            }
            self.assembly['instructions'].append(op_instruction)
            self.instruction_id += 1
            return

        if gate_key in self.symbol_table['gates']:
            implementation = self.symbol_table['gates'][gate_key]
            self.visit(implementation)
        else:
            self.error_reporter.report(f"Gate '{gate_name}' on {targets} is not defined.", node.line, node.column)

    def visit_MeasureOperation(self, node):
        qubit_name = node.qubit.name
        if qubit_name not in self.symbol_table['qubits']:
            self.error_reporter.report(f"Qubit '{qubit_name}' is not defined.", node.line, node.column)
            return

        classical_bit_name = node.classical_bit.name

        measure_instruction = {
            "instruction_id": self.instruction_id,
            "command": "MEASURE",
            "params": {
                "qubit": qubit_name,
                "classical_bit": classical_bit_name
            }
        }
        self.assembly['instructions'].append(measure_instruction)
        self.instruction_id += 1

    def visit_MeasureStatement(self, node):
        qubit_name = node.qubit.name
        classical_bit_name = node.classical_bit.name if node.classical_bit else None
        
        if qubit_name not in self.symbol_table['qubits']:
            self.error_reporter.report(f"Qubit '{qubit_name}' used in measure statement is not defined in this scope.", node.line, node.column)
        
        if classical_bit_name and classical_bit_name not in self.symbol_table['classical_bits']:
            self.error_reporter.report(f"Classical bit '{classical_bit_name}' used in measure statement is not defined in this scope.", node.line, node.column)

        measure_instruction = {
            "instruction_id": self.instruction_id,
            "command": "MEASURE",
            "params": {
                "qubit": qubit_name,
                "classical_bit": classical_bit_name
            }
        }
        self.assembly['instructions'].append(measure_instruction)
        self.instruction_id += 1

    def visit_ClassicalBitDefinition(self, node):
        bit_name = node.name.name
        if bit_name in self.symbol_table['classical_bits']:
            self.error_reporter.report(f"Classical bit '{bit_name}' is already defined in this scope.", node.line, node.column)
        self.symbol_table['classical_bits'][bit_name] = {'line': node.line, 'column': node.column}

    def visit_QubitDeclaration(self, node):
        name = node.name.name
        if name in self.symbol_table['qubits']:
            self.error_reporter.report(f"Qubit '{name}' is already defined in this scope.", node.line, node.column)
        self.symbol_table['qubits'][name] = {'line': node.line, 'column': node.column}

def main():
    if len(sys.argv) != 2:
        print("Usage: python -m arbol.compiler.compiler <path_to_arbol_file>")
        sys.exit(1)

    arbol_file_path = sys.argv[1]
    output_json_path = arbol_file_path.replace('.arbol', '.bsim.json')

    try:
        with open(arbol_file_path, 'r') as f:
            arbol_code = f.read()
        
        error_reporter = ErrorReporter()
        lexer = Lexer(arbol_code)
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()

        if error_reporter.has_errors():
            error_reporter.display()
            sys.exit(1)

        compiler = Compiler(error_reporter)
        assembly = compiler.compile(ast)

        if error_reporter.has_errors():
            error_reporter.display()
            sys.exit(1)

        with open(output_json_path, 'w') as f:
            json.dump(assembly, f, indent=4)

        print(f"Successfully compiled {arbol_file_path} to {output_json_path}")

    except FileNotFoundError:
        print(f"Error: File not found at {arbol_file_path}")
    except Exception as e:
        print(f"An unexpected error occurred during compilation: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
