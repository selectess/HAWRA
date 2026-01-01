import sys
import os
import json
from lark import Lark, Transformer, v_args
from .arbol_ast import *
from .error import ErrorReporter, CompilationError

# Grammaire Arbol unifiée utilisant Lark
ARBOL_GRAMMAR = r"""
    ?start: program

    program: (statement ";"?)*

    ?statement: gate_definition
              | stimulus_definition
              | circuit_definition
              | config_definition
              | run_statement
              | genes_block
              | grn_definition
              | logical_qubit_definition
              | qubit_declaration
              | gate_application
              | measure_statement

    qubit_declaration: "qubit" CNAME
    gate_application: CNAME CNAME ("," CNAME)*
    measure_statement: "measure" CNAME ["->" CNAME]

    gate_definition: "gate" CNAME "(" [param_list] ")" "{" quantum_operation* "}"
    param_list: CNAME ("," CNAME)*

    stimulus_definition: "stimulus" CNAME "(" [param_def_list] ")" "{" "}"
    param_def_list: param_def ("," param_def)*
    param_def: CNAME ":" CNAME

    circuit_definition: "circuit" CNAME "(" [param_def_list] ")" "{" circuit_body_statement* "}"
    ?circuit_body_statement: (qubit_def | bit_def | quantum_operation | stimulus_application | measure_assignment) ";"?
    qubit_def: CNAME ":" "qubit"
    bit_def: CNAME ":" "bit"

    quantum_operation: CNAME "(" [param_list] ")"

    stimulus_application: "apply" CNAME "(" [arg_list] ")" "on" CNAME
    arg_list: arg ("," arg)*
    arg: CNAME ":" value
    ?value: SIGNED_NUMBER | STRING | CNAME

    measure_assignment: CNAME "=" "measure" CNAME

    config_definition: "config" "{" [config_section ("," config_section)*] "}"
    config_section: CNAME "{" [config_entry ("," config_entry)*] ","? "}"
    config_entry: CNAME ":" value

    run_statement: "run" CNAME "(" [arg_list] ")"
                 | "step" value

    genes_block: "genes" "{" gene_definition* "}"
    gene_definition: CNAME "{" config_entry* "}"

    grn_definition: "grn" "{" grn_interaction* "}"
    grn_interaction: CNAME ("activates" | "represses") CNAME ("[" arg_list "]")?

    logical_qubit_definition: "LOGICAL_QUBIT" CNAME ("[" SIGNED_NUMBER "]")? "is" "{" promoter_list "}"
    promoter_list: promoter ("," promoter)*
    promoter: CNAME ":" CNAME

    COMMENT_MULTI: /\(\*[\s\S]*?\*\)/
    COMMENT_SINGLE: /\/\/[^\n]*/
    %ignore COMMENT_MULTI
    %ignore COMMENT_SINGLE

    %import common.CNAME
    %import common.ESCAPED_STRING -> STRING
    %import common.SIGNED_NUMBER
    %import common.WS
    %ignore WS
"""

class ArbolTransformer(Transformer):
    def __init__(self, compiler):
        super().__init__()
        self.compiler = compiler

    def CNAME(self, cname): return str(cname)
    def STRING(self, s): return s[1:-1]
    def SIGNED_NUMBER(self, n): return float(n)

    def program(self, statements): return statements
    
    def circuit_definition(self, items):
        name = items[0]
        # param_def_list est une liste de Parameter
        params = items[1] if items[1] is not None else []
        body = items[2:]
        return CircuitDefinition(name=Identifier(name), parameters=params, body=body)

    def param_def_list(self, items):
        return items

    def param_def(self, items):
        return Parameter(name=items[0], type=items[1])

    def qubit_def(self, items):
        return QubitDeclaration(name=Identifier(items[0]))

    def bit_def(self, items):
        # On peut réutiliser ClassicalBitDefinition ou créer un BitDeclaration
        return ClassicalBitDefinition(name=Identifier(items[0]), line=0, column=0)

    def quantum_operation(self, items):
        name = items[0]
        # param_list est une liste de CNAME (noms de qubits ou params)
        refs = [Identifier(q) for q in items[1]] if len(items) > 1 and items[1] else []
        return GateApplication(gate_name=Identifier(name), arguments=[], qubits=refs, line=0, column=0)

    def param_list(self, items):
        return items

    def stimulus_application(self, items):
        name = items[0]
        args = items[1]
        target = items[2]
        return StimulusApplication(stimulus_name=Identifier(name), arguments=args, target=Identifier(target))

    def measure_assignment(self, items):
        target_bit = items[0]
        qubit = items[1]
        return MeasureStatement(qubit=Identifier(qubit), classical_bit=Identifier(target_bit), line=0, column=0)

    def config_entry(self, items):
        name, val = items
        return ConfigParameter(name=Identifier(name), value=val)

    def gene_definition(self, items):
        name = items[0]
        props = items[1:]
        return GeneDefinition(gene_name=Identifier(name), properties=props)

    def genes_block(self, items):
        return GenesBlock(genes=items)

    def grn_interaction(self, items):
        reg, type_str, target = items[0:3]
        params = items[3] if len(items) > 3 else {}
        param_nodes = [ConfigParameter(Identifier(k), v) for k, v in params.items()]
        return GRNInteraction(
            regulator=Identifier(reg),
            interaction_type="ACTIVATES" if "activates" in type_str else "REPRESSES",
            target=Identifier(target),
            params=param_nodes
        )

    def grn_definition(self, items):
        return GRNDefinition(interactions=items)

    def run_statement(self, items):
        if items[0] == "step":
            return StepInstruction(duration=Argument(items[1]))
        return RunCircuit(circuit_name=Identifier(items[0]), arguments=items[1])

    def arg_list(self, items):
        return {k: v for k, v in items}

    def arg(self, items):
        return items[0], items[1]

    def qubit_declaration(self, items):
        return QubitDeclaration(name=Identifier(items[0]))

    def gate_application(self, items):
        gate_name = items[0]
        qubits = [Identifier(q) for q in items[1:]]
        return GateApplication(gate_name=Identifier(gate_name), arguments=[], qubits=qubits, line=0, column=0)

    def measure_statement(self, items):
        qubit = items[0]
        cbit = items[1] if len(items) > 1 else None
        return MeasureStatement(qubit=Identifier(qubit), classical_bit=Identifier(cbit) if cbit else None, line=0, column=0)

class Compiler:
    def __init__(self, error_reporter=None):
        self.assembly = {
            "metadata": {
                "source_arbol": "compiled.arbol",
                "version": "0.3",
                "qubit_count": 0,
                "qubits": []
            },
            "instructions": []
        }
        self.instruction_id = 0
        self.symbol_table = {'genes': {}, 'stimuli': {}, 'circuits': {}, 'qubits': {}}
        self.error_reporter = error_reporter if error_reporter else ErrorReporter()
        self.initialize_assembly()

    def initialize_assembly(self):
        self.add_instruction("INITIALIZE", {
            "max_time": 500,
            "dt": 1.0,
            "env": {"light_schedule": []},
            "quantum": {"p700_threshold": 0.8, "decoherence_rate": 0.015},
            "bio": {
                "p700_synthesis_rate": 0.05,
                "p700_degradation_rate": 0.1,
                "optimal_temp": 35.0,
                "genes": {},
                "grn": {}
            }
        })

    def add_instruction(self, command, params):
        self.assembly['instructions'].append({
            "instruction_id": self.instruction_id,
            "command": command,
            "params": params if command != "INITIALIZE" else None,
            "config": params if command == "INITIALIZE" else None
        })
        self.instruction_id += 1

    def compile(self, source):
        if isinstance(source, str):
            from .lexer import Lexer
            from .parser import Parser
            lexer = Lexer(source)
            parser = Parser(lexer, self.error_reporter)
            ast = parser.parse()
            if self.error_reporter.has_errors():
                return None
            return self.compile(ast)
        
        # Si c'est déjà un AST
        try:
            if hasattr(source, 'statements'):
                for stmt in source.statements:
                    self.visit(stmt)
            else:
                self.visit(source)
            return self.assembly
        except Exception as e:
            self.error_reporter.report(f"Compilation Error: {str(e)}", 0, 0)
            return None

    def visit(self, node):
        method_name = 'visit_' + type(node).__name__
        visitor = getattr(self, method_name, self.generic_visit)
        return visitor(node)

    def generic_visit(self, node):
        pass

    def visit_GenesBlock(self, node):
        for gene_def in node.genes:
            name = gene_def.gene_name.name
            props = {p.name.name: p.value for p in gene_def.properties}
            self.assembly['instructions'][0]['config']['bio']['genes'][name] = props

    def visit_GRNDefinition(self, node):
        grn = self.assembly['instructions'][0]['config']['bio']['grn']
        for inter in node.interactions:
            target = inter.target.name
            if target not in grn: grn[target] = {"regulators": []}
            
            params = {p.name.name: p.value for p in inter.params}
            grn[target]["regulators"].append({
                "gene": inter.regulator.name,
                "type": "activator" if inter.interaction_type == "ACTIVATES" else "repressor",
                "weight": params.get('weight', 1.0),
                "hill_coefficient": params.get('hill_coefficient', 2.0),
                "half_max_concentration": params.get('half_max_concentration', 0.5)
            })

    def visit_StepInstruction(self, node):
        self.add_instruction("RUN_UNTIL", {"time": node.duration.value})

    def visit_CircuitDefinition(self, node):
        self.symbol_table['circuits'][node.name.name] = node

    def visit_RunCircuit(self, node):
        name = node.circuit_name.name
        if name in self.symbol_table['circuits']:
            # Marqueur de début de run pour le BioOS et les métriques
            self.add_instruction("RUN_START", {"circuit": name, "args": node.arguments})
            
            circuit_node = self.symbol_table['circuits'][name]
            # Expansion du corps du circuit
            for stmt in circuit_node.body:
                self.visit(stmt)
                
            self.add_instruction("RUN_END", {"circuit": name})
        else:
            self.error_reporter.report(f"Circuit '{name}' non défini.", 0, 0)

    def visit_StimulusApplication(self, node):
        self.add_instruction("STIMULUS_APPLY", {
            "stimulus": node.stimulus_name.name,
            "arguments": node.arguments,
            "target": node.target.name
        })

    def visit_QubitDeclaration(self, node):
        name = node.name.name
        if name not in self.symbol_table['qubits']:
            self.assembly['metadata']['qubit_count'] += 1
            self.assembly['metadata']['qubits'].append({"name": name})
            self.symbol_table['qubits'][name] = self.assembly['metadata']['qubit_count'] - 1

    def visit_GateApplication(self, node):
        gate = node.gate_name.name.lower()
        targets = []
        for q_id in node.qubits:
            q_name = q_id.name
            if q_name not in self.symbol_table['qubits']:
                # Auto-déclaration si manquant pour les tests rapides
                self.visit_QubitDeclaration(QubitDeclaration(name=q_id))
            
            # Utiliser l'index au lieu du nom pour le simulateur
            targets.append(self.symbol_table['qubits'][q_name])
        
        self.add_instruction("GATE_APPLY", {"gate": gate, "targets": targets})

    def visit_MeasureStatement(self, node):
        qubit = node.qubit.name
        if qubit not in self.symbol_table['qubits']:
            self.visit_QubitDeclaration(QubitDeclaration(name=node.qubit))
        
        target_idx = self.symbol_table['qubits'][qubit]
        cbit = node.classical_bit.name if node.classical_bit else "c0"
        self.add_instruction("MEASURE", {"qubit": target_idx, "classical_bit": cbit})

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python -m arbol.compiler.compiler <source.arbol>")
        sys.exit(1)
    
    with open(sys.argv[1], 'r') as f:
        source = f.read()
    
    reporter = ErrorReporter()
    compiler = Compiler(reporter)
    result = compiler.compile(source)
    
    if result:
        output_file = sys.argv[1].replace('.arbol', '.bsim.json')
        with open(output_file, 'w') as f:
            json.dump(result, f, indent=4)
    else:
        reporter.display()
        sys.exit(1)
