import json
from lark import Lark, Transformer, v_args

arbol_grammar = r"""
    ?start: program

    program: (statement ";"?)*

    ?statement: gate_definition
              | stimulus_definition
              | circuit_definition
              | config_definition
              | run_statement

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

class ArbolParser:
    def __init__(self):
        self.parser = Lark(arbol_grammar, start='program', parser='lalr')

    def parse(self, text):
        return self.parser.parse(text)

@v_args(inline=True)
class ArbolTransformer(Transformer):
    def __init__(self):
        super().__init__()
        self.instructions = []
        self.instruction_id = 0
        self.circuits = {}
        self.qubits = set()
        self.bits = set()

    def program(self, *statements):
        # The run_statement will populate the instructions
        # We only care about the side-effects of the transformation
        pass

    def CNAME(self, cname):
        return cname.value

    def STRING(self, s):
        return s[1:-1].replace('\\"', '"')

    def SIGNED_NUMBER(self, n):
        return float(n)

    def value(self, v):
        return v

    def param_list(self, *params):
        return list(params)

    def arg(self, name, value):
        return (name, value)

    def arg_list(self, *args):
        return dict(args)

    def param_def(self, name, type):
        return (name, type)

    def param_def_list(self, *params):
        return dict(params)

    def gate_definition(self, name, params, *body):
        # For now, we don't need to handle gate definitions for the simulation
        pass

    def stimulus_definition(self, name, params):
        # For now, we don't need to handle stimulus definitions for the simulation
        pass

    def circuit_definition(self, name, params, *body):
        self.circuits[name] = {'params': params or {}, 'body': body}

    def qubit_def(self, name):
        self.qubits.add(name)
        return ("qubit_def", name)

    def bit_def(self, name):
        self.bits.add(name)
        return ("bit_def", name)

    def quantum_operation(self, name, params):
        return ("quantum_operation", name, params or [])

    def stimulus_application(self, name, args, target):
        return ("stimulus_application", name, args or {}, target)

    def measure_assignment(self, cbit, qbit):
        return ("measure", qbit, cbit)

    def circuit_body_statement(self, stmt):
        return stmt

    def run_statement(self, circuit_name, args):
        if circuit_name not in self.circuits:
            raise Exception(f"Circuit '{circuit_name}' not defined.")

        circuit = self.circuits[circuit_name]
        # We can handle arguments later if needed.

        for stmt in circuit['body']:
            if not stmt:
                continue
            
            if isinstance(stmt, tuple):
                stmt_type = stmt[0]
                if stmt_type == "quantum_operation":
                    _, gate_name, target_qubits = stmt
                    instruction = {
                        "op_code": "APPLY_GATE",
                        "params": {
                            "gate": gate_name,
                            "target": target_qubits
                        }
                    }
                    self.instructions.append(instruction)
                elif stmt_type == "measure":
                    _, source_qubit, target_bit = stmt
                    instruction = {
                        "op_code": "MEASURE",
                        "params": {
                            "target": [source_qubit],
                            "cbit": [target_bit]
                        }
                    }
                    self.instructions.append(instruction)
                elif stmt_type == "qubit_def":
                    _, name = stmt
                    self.qubits.add(name)
    
    def string(self, s):
        return s[1:-1].replace('\\"', '"')

class ArbolCompiler:
    def __init__(self, source_code):
        self.source_code = source_code
        self.tree = ArbolParser().parse(source_code)
        self.transformer = ArbolTransformer()
        self.transformer.transform(self.tree)
        self.instructions = self.transformer.instructions
        self.qubits = list(self.transformer.qubits)

    def to_json(self, source_filename="compiled.arbol"):
        output = {
            "metadata": {
                "source_arbol": source_filename,
                "version": "0.2",
                "qubit_count": len(self.qubits),
                "qubits": self.qubits
            },
            "instructions": self.instructions
        }
        return json.dumps(output, indent=2)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        source_path = sys.argv[1]
        with open(source_path, 'r') as f:
            source_code = f.read()
        compiler = ArbolCompiler(source_code)
        output_path = source_path.replace('.arbol', '.bsim.json')
        with open(output_path, 'w') as f:
            f.write(compiler.to_json(source_filename=source_path.split('/')[-1]))
        print(f"Compiled {source_path} to {output_path}")
    else:
        print("Usage: python arbol_compiler.py <source_file.arbol>")