from .lexer import Lexer, Token
from .arbol_ast import (
    Program, CircuitDefinition, GateDefinition, QuantumOperation, MeasureOperation,
    QubitDeclaration, GateApplication, StimulusDefinition, LogicalQubitDefinition,
    StimulusApplication, MeasureStatement, ClassicalBitDefinition, ConfigDefinition, ConfigSection,
    ConfigParameter, RunCircuit, GeneDefinition, GenesBlock, GRNInteraction,
    GRNDefinition, StepInstruction, RunBlock, Identifier, Parameter, Argument,
    Promoter, IsClause
)
from .error import CompilationError

class Parser:
    def __init__(self, lexer, error_reporter):
        self.lexer = lexer
        self.tokens = list(lexer.tokens())  # Convertir le générateur en liste
        self.pos = 0
        self.error_reporter = error_reporter

    def parse(self):
        print("[DEBUG] --- Calling Parser.parse() ---")
        statements = []
        while self.current_token is not None:
            try:
                statement = self.parse_statement()
                if statement:
                    statements.append(statement)
                else:
                    break
            except CompilationError as e:
                self.error_reporter.report(e.message, e.line, e.column)
                self.synchronize()
        return Program(statements=statements)

    @property
    def current_token(self):
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def advance(self):
        self.pos += 1

    def peek_token(self):
        next_pos = self.pos + 1
        if next_pos < len(self.tokens):
            return self.tokens[next_pos]
        return None

    def eat(self, token_type):
        token = self.current_token
        if token is None:
            raise CompilationError(f"Expected {token_type}, but found end of file.", 1, 1)
        if token.type != token_type:
            raise CompilationError(f"Expected {token_type} but got {token.type}", token.line, token.column)
        self.advance()
        return token

    def synchronize(self):
        while self.current_token:
            if self.current_token.type == 'SEMICOLON':
                self.advance()
                return
            if self.current_token.type in ['CIRCUIT', 'GATE', 'STIMULUS', 'CONFIG', 'LOGICAL_QUBIT', 'RUN']:
                return
            self.advance()

    def parse_statement(self):
        token = self.current_token
        if not token:
            return None

        if token.type == 'CIRCUIT':
            return self.parse_circuit_definition()
        elif token.type == 'GATE':
            return self.parse_gate_definition()
        elif token.type == 'LOGICAL_QUBIT':
            print(f"[DEBUG] Found LOGICAL_QUBIT token: {token}")
            return self.parse_logical_qubit_definition()
        elif token.type == 'STIMULUS':
            return self.parse_stimulus_definition()
        elif token.type == 'CLASSICAL_BIT':
            return self.parse_classical_bit_definition()
        elif token.type == 'MEASURE':
            return self.parse_measure_statement()
        elif token.type == 'APPLY':
            return self.parse_stimulus_application()
        elif token.type == 'GENES':
            return self.parse_genes_definition()
        elif token.type == 'GRN':
            return self.parse_grn_definition()
        elif token.type == 'RUN':
            return self.parse_run_block()
        elif token.type == 'CONFIG':
            return self.parse_config_definition()
        elif token.type in ['H', 'X', 'Y', 'Z', 'CNOT', 'CCNOT', 'DJ_ORACLE_CONSTANT', 'DJ_ORACLE_BALANCED', 'ID', 'CX']:
            return self.parse_quantum_gate()
        else:
            if token and hasattr(token, 'type') and hasattr(token, 'line') and hasattr(token, 'column'):
                raise CompilationError(f"Unexpected token {token.type}", token.line, token.column)
            else:
                raise CompilationError("Unexpected token (unknown type)", 1, 1)

    def parse_quantum_gate(self):
        gate_token = self.current_token
        if gate_token is None:
            raise CompilationError("Expected quantum gate but found end of file", 1, 1)
        gate_name = self.eat(gate_token.type).value
        qubits = []
        
        # Handle H(q1) syntax
        if self.current_token and self.current_token.type == 'LPAREN':
            self.eat('LPAREN')
            qubits = self.parse_identifier_list()
            self.eat('RPAREN')
        # Handle H on q1 syntax
        elif self.current_token and self.current_token.type == 'ON':
            self.eat('ON')
            while self.current_token and self.current_token.type == 'ID':
                qubits.append(Identifier(self.eat('ID').value))
        # Handle H q1 syntax
        elif self.current_token and self.current_token.type == 'ID':
             while self.current_token and self.current_token.type == 'ID':
                qubits.append(Identifier(self.eat('ID').value))
        
        self.eat('SEMICOLON')
        return GateApplication(gate_name=Identifier(gate_name), arguments=[], qubits=qubits, line=gate_token.line, column=gate_token.column)



























    def parse_config_definition(self):
        self.eat('CONFIG')
        self.eat('LBRACE')
        sections = []
        while self.current_token and self.current_token.type != 'RBRACE':
            section_name = self.eat('ID').value
            self.eat('LBRACE')
            params = []
            while self.current_token and self.current_token.type != 'RBRACE':
                param_name = self.eat('ID').value
                self.eat('COLON')
                if self.current_token.type == 'NUMBER':
                    param_val = float(self.eat('NUMBER').value)
                elif self.current_token.type == 'STRING':
                    param_val = self.eat('STRING').value
                else:
                    # Fallback or error - assume identifier or keyword
                    token = self.current_token
                    param_val = token.value
                    self.advance() # consume whatever it is
                
                params.append(ConfigParameter(name=Identifier(name=param_name), value=param_val))
                if self.current_token.type == 'COMMA':
                    self.eat('COMMA')
            self.eat('RBRACE')
            sections.append(ConfigSection(name=Identifier(name=section_name), parameters=params))
            if self.current_token and self.current_token.type == 'COMMA':
                self.eat('COMMA')
        self.eat('RBRACE')
        return ConfigDefinition(sections=sections)

    def parse_logical_qubit_definition(self):
        try:
            self.eat('LOGICAL_QUBIT')
            name_token = self.eat('ID')
            name = Identifier(name_token.value)
            
            if self.current_token and self.current_token.type == 'LBRACKET':
                self.eat('LBRACKET')
                size = int(self.eat('NUMBER').value)
                self.eat('RBRACKET')
            else:
                size = 1
            
            self.eat('IS')
            
            promoters = []
            self.eat('LBRACE')
            while self.current_token and self.current_token.type != 'RBRACE':
                promoter_type_token = self.current_token
                if promoter_type_token.type not in ['ID', 'ACTIVATOR', 'REPRESSOR']:
                    raise CompilationError(f"Expected promoter type (activator/repressor) but got {promoter_type_token.type}", promoter_type_token.line, promoter_type_token.column)
                self.advance()
                promoter_type = promoter_type_token.value

                if promoter_type.lower() not in ['activator', 'repressor']:
                    raise CompilationError(f"Unknown promoter type '{promoter_type}'", promoter_type_token.line, promoter_type_token.column)
                
                self.eat('COLON')
                promoter_name_token = self.eat('ID')
                promoter_name = promoter_name_token.value
                promoters.append(Promoter(type=promoter_type, name=promoter_name))
                
                if self.current_token and self.current_token.type == 'COMMA':
                    self.eat('COMMA')
                elif self.current_token and self.current_token.type != 'RBRACE':
                    raise CompilationError(f"Expected COMMA or RBRACE but got {self.current_token.type}", self.current_token.line, self.current_token.column)

            self.eat('RBRACE')
            self.eat('SEMICOLON')
            
            is_clause = IsClause(promoters=promoters)
            result = LogicalQubitDefinition(name=name, size=size, is_clause=is_clause, line=name_token.line, column=name_token.column)
            return result
        except CompilationError as e:
            raise

    def parse_stimulus_definition(self):
        self.eat('STIMULUS')
        name = Identifier(self.eat('ID').value)
        self.eat('LPAREN')
        parameters = []
        while self.current_token and self.current_token.type != 'RPAREN':
            param_name = self.eat('ID').value
            self.eat('COLON')
            if self.current_token and self.current_token.type == 'STRING':
                param_default = self.eat('STRING').value
                param_type = 'STRING'
            elif self.current_token and self.current_token.type == 'NUMBER':
                param_default = float(self.eat('NUMBER').value)
                param_type = 'NUMBER'
            elif self.current_token and self.current_token.type == 'ID':
                # Handle type specifier (e.g. intensity: float, type: light)
                param_default = self.eat('ID').value
                param_type = param_default # Treat the ID as the type/default
            else:
                param_default = None
                param_type = 'STRING'
            parameters.append(Parameter(name=param_name, type=param_type))
            if self.current_token and self.current_token.type == 'COMMA':
                self.eat('COMMA')
            else:
                break
        self.eat('RPAREN')
        
        # Optional body { ... }
        if self.current_token and self.current_token.type == 'LBRACE':
            self.eat('LBRACE')
            # consume body until RBRACE (naively for now, just counting braces or skipping)
            # Since we don't have a parse_block that skips, we'll just skip tokens until RBRACE
            brace_count = 1
            while brace_count > 0 and self.current_token:
                if self.current_token.type == 'LBRACE':
                    brace_count += 1
                elif self.current_token.type == 'RBRACE':
                    brace_count -= 1
                self.advance()
            # The last advance consumed the matching RBRACE
        else:
            self.eat('SEMICOLON')
            
        return StimulusDefinition(name=name, parameters=parameters, promoter=Promoter(type='dummy', name='dummy'), is_clause=IsClause(promoters=[]), line=0, column=0)

    def parse_classical_bit_definition(self):
        token = self.current_token
        if not token:
            raise CompilationError("Unexpected end of file while parsing classical bit definition.", 1, 1)
        self.eat('CLASSICAL_BIT')
        name = Identifier(self.eat('ID').value)
        self.eat('SEMICOLON')
        return ClassicalBitDefinition(name=name, line=token.line, column=token.column)



    def parse_is_clause(self):
        self.eat('LBRACE')
        promoters = []
        while self.current_token and self.current_token.type != 'RBRACE':
            promoter_type = self.current_token.value
            self.advance()
            self.eat('COLON')
            promoter_name = self.eat('ID').value
            promoters.append(Promoter(type=promoter_type, name=promoter_name))
            if self.current_token.type == 'RBRACE':
                break
            self.eat('COMMA')
        self.eat('RBRACE')
        return IsClause(promoters=promoters)

    def parse_genes_definition(self):
        self.eat('GENES')
        self.eat('LBRACE')
        genes = []
        while self.current_token and self.current_token.type != 'RBRACE':
            self.eat('ID')  # gene
            gene_name = self.eat('STRING').value
            properties = []
            if self.current_token.type == 'WITH':
                self.eat('WITH')
                while self.current_token.type != 'SEMICOLON':
                    prop_name = self.eat(self.current_token.type).value
                    self.eat('EQUALS')
                    prop_value = self.eat('NUMBER').value
                    try:
                        prop_value = float(prop_value)
                    except ValueError:
                        pass  # Keep as string if conversion fails
                    properties.append(ConfigParameter(name=Identifier(name=prop_name), value=prop_value))
                    if self.current_token and self.current_token.type == 'COMMA':
                        self.eat('COMMA')
                    else:
                        break  # finished with properties
            genes.append(GeneDefinition(gene_name=Identifier(name=gene_name), properties=properties))
            self.eat('SEMICOLON')
        self.eat('RBRACE')
        return GenesBlock(genes=genes)

    def parse_grn_definition(self):
        self.eat('GRN')
        self.eat('LBRACE')
        interactions = []
        while self.current_token and self.current_token.type != 'RBRACE':
            source_gene = self.eat('STRING').value
            if self.current_token.type in ['ACTIVATES', 'REPRESSES']:
                interaction_type = self.eat(self.current_token.type).type
            else:
                self.error_reporter.report(f"Unexpected token {self.current_token.type}", self.current_token.line, self.current_token.column)
                interaction_type = '' # Initialize in case of error
            target_gene = self.eat('STRING').value
            properties = []
            if self.current_token.type == 'WITH':
                self.eat('WITH')
                while self.current_token.type != 'SEMICOLON':
                    prop_name = self.eat(self.current_token.type).value
                    self.eat('EQUALS')
                    prop_value = self.eat('NUMBER').value
                    properties.append(ConfigParameter(name=Identifier(name=prop_name), value=prop_value))
                    if self.current_token and self.current_token.type == 'COMMA':
                        self.eat('COMMA')
                    else:
                        break
            interactions.append(GRNInteraction(
                regulator=Identifier(name=source_gene),
                target=Identifier(name=target_gene),
                interaction_type=interaction_type,
                params=properties
            ))
            self.eat('SEMICOLON')
        self.eat('RBRACE')
        return GRNDefinition(interactions=interactions)

    def parse_run_block(self):
        run_token = self.current_token
        self.eat('RUN')
        if self.current_token and self.current_token.type == 'LBRACE':
            self.eat('LBRACE')
            instructions = []
            while self.current_token and self.current_token.type != 'RBRACE':
                if self.current_token.type == 'STEP':
                    self.eat('STEP')
                    self.eat('LPAREN')
                    duration_arg = Argument(value=self.eat('NUMBER').value)
                    self.eat('RPAREN')
                    instructions.append(StepInstruction(duration=duration_arg))
                elif self.current_token and self.current_token.type == 'APPLY':
                    instructions.append(self.parse_stimulus_application())
                self.eat('SEMICOLON')
            self.eat('RBRACE')
            return RunBlock(instructions=instructions)
        else:
            # Parse run circuit_name(args);
            circuit_name = Identifier(self.eat('ID').value)
            arguments = []
            if self.current_token and self.current_token.type == 'LPAREN':
                self.eat('LPAREN')
                while self.current_token and self.current_token.type != 'RPAREN':
                    arg_name = self.eat('ID').value
                    self.eat('COLON')
                    arg_value = self.eat('NUMBER').value # Assuming numbers for now
                    arguments.append(ConfigParameter(name=Identifier(name=arg_name), value=float(arg_value)))
                    if self.current_token.type == 'COMMA':
                        self.eat('COMMA')
                    else:
                        break
                self.eat('RPAREN')
            self.eat('SEMICOLON')
            return RunCircuit(circuit_name=circuit_name, arguments=arguments, line=run_token.line if run_token else 0, column=run_token.column if run_token else 0)

    def parse_stimulus_application(self):
        self.eat('APPLY')
        stimulus_name = self.eat('ID').value
        self.eat('LPAREN')
        args = []
        while self.current_token and self.current_token.type != 'RPAREN':
            arg_name = self.eat('ID').value
            self.eat('COLON')
            arg_value = self.eat('NUMBER').value
            args.append(ConfigParameter(name=Identifier(name=arg_name), value=arg_value))
            if self.current_token.type == 'COMMA':
                self.eat('COMMA')
            else:
                break
        self.eat('RPAREN')
        self.eat('ON')
        if self.current_token and self.current_token.type == 'STRING':
            target_gene = self.eat('STRING').value
        else:
            target_gene = self.eat('ID').value
        return StimulusApplication(stimulus_name=Identifier(name=stimulus_name), arguments=args, target=Identifier(name=target_gene))

    def parse_circuit_definition(self):
        self.eat('CIRCUIT')
        name = Identifier(self.eat('ID').value)
        parameters = self.parse_parameters()
        self.eat('LBRACE')
        body = []
        while self.current_token and self.current_token.type != 'RBRACE':
            statement = self.parse_circuit_statement()
            if statement:
                body.append(statement)
        self.eat('RBRACE')
        return CircuitDefinition(name=name, parameters=parameters, body=body)

    def parse_gate_definition(self):
        self.eat('GATE')
        token = self.current_token
        if not token:
            raise CompilationError("Unexpected end of file after GATE keyword.", 1, 1)
        if token.type == 'ID' or token.type in ['H', 'X', 'Y', 'Z', 'CNOT']:
            name = Identifier(self.eat(token.type).value)
        else:
            raise CompilationError(f"Expected gate name (Identifier) but got {token.type}", token.line, token.column)
        
        # Handle both gate H(q1) and gate H on q1 syntax
        if self.current_token and self.current_token.type == 'LPAREN':
            self.eat('LPAREN')
            qubit_refs = self.parse_identifier_list()
            self.eat('RPAREN')
            self.eat('LBRACE')
            body = []
            while self.current_token and self.current_token.type != 'RBRACE':
                op = self.parse_quantum_operation()
                if op:
                    body.append(op)
            self.eat('RBRACE')
            return GateDefinition(name=name, qubit_refs=qubit_refs, body=body)
        elif self.current_token and self.current_token.type == 'ON':
            self.eat('ON')
            # For single qubit case like "on q1 is", parse manually
            qubit_id = self.eat('ID').value
            qubit_refs = [Identifier(qubit_id)]
            self.eat('IS')
            # Parse the 'is' clause - handle the stimulus application
            body = []
            if self.current_token and self.current_token.type == 'APPLY':
                # Parse apply light to q1 with { ... };
                self.eat('APPLY')
                stimulus_name = self.eat('ID').value
                self.eat('TO')
                target_qubit = self.eat('ID').value
                self.eat('WITH')
                self.eat('LBRACE')
                
                # Parse the parameters
                params = {}
                while self.current_token and self.current_token.type != 'RBRACE':
                    param_name = self.eat('ID').value
                    self.eat('COLON')
                    param_value = self.eat('STRING').value
                    params[param_name] = param_value
                    if self.current_token and self.current_token.type == 'COMMA':
                        self.eat('COMMA')
                    else:
                        break
                self.eat('RBRACE')
                self.eat('SEMICOLON')
                
                stim_app = StimulusApplication(
                    stimulus_name=Identifier(name=stimulus_name),
                    arguments=[ConfigParameter(name=Identifier(k), value=v) for k, v in params.items()],
                    target=Identifier(name=target_qubit)
                )
                body.append(stim_app)
            return GateDefinition(name=name, qubit_refs=qubit_refs, body=body)
        else:
            # Fallback or error for when neither ( nor ON is present
            token = self.current_token
            if token:
                raise CompilationError(f"Expected '(' or 'on' after gate name, but got {token.type}", token.line, token.column)
            else:
                raise CompilationError("Unexpected end of file after gate name.", 1, 1)


    def parse_parameters(self):
        parameters = []
        if self.current_token and self.current_token.type == 'LPAREN':
            self.eat('LPAREN')
            if self.current_token and self.current_token.type == 'RPAREN':
                self.eat('RPAREN')
                return parameters

            while self.current_token and self.current_token.type != 'RPAREN':
                name = self.eat('ID').value
                self.eat('COLON')
                param_type = self.eat('ID').value
                
                if param_type == 'variable':
                    self.eat('LPAREN')
                    self.eat('NUMBER')
                    self.eat('TO')
                    self.eat('NUMBER')
                    self.eat('OVER')
                    self.eat('NUMBER')
                    self.eat('RPAREN')

                parameters.append(Parameter(name=name, type=param_type))
                
                if self.current_token and self.current_token.type == 'COMMA':
                    self.eat('COMMA')
                elif self.current_token and self.current_token.type != 'RPAREN':
                    token = self.current_token
                    raise CompilationError(f"Expected COMMA or RPAREN but got {token.type}", token.line, token.column)

            self.eat('RPAREN')
        return parameters

    def parse_identifier_list(self):
        identifiers = []
        stop_tokens = ['RPAREN', 'RBRACE']
        
        # Check if current_token is None right at the start
        if not self.current_token:
            return identifiers

        if self.current_token.type not in stop_tokens:
            while self.current_token and self.current_token.type not in stop_tokens:
                identifiers.append(Identifier(self.eat('ID').value))
                if self.current_token and self.current_token.type == 'COMMA':
                    self.eat('COMMA')
                else:
                    break  # Exit if no comma or if token is None
        return identifiers

    def parse_qubit_list(self):
        qubits = []
        # Ensure there is a token and it's an ID
        while self.current_token and self.current_token.type == 'ID':
            qubits.append(Identifier(self.eat('ID').value))
            # If there's a comma, consume it and expect another qubit
            if self.current_token and self.current_token.type == 'COMMA':
                self.eat('COMMA')
            else:
                break  # Exit loop if there is no comma
        return qubits

    def parse_circuit_statement(self):
        token = self.current_token
        if not token:
            return None  # End of file

        peek = self.peek_token()
        if not peek:
            raise CompilationError(f"Unexpected end of file after {token.type}", token.line, token.column)

        if token.type == 'ID' and peek.type == 'COLON':
            return self.parse_qubit_declaration()
        elif token.type == 'ID' and peek.type == 'LPAREN':
            return self.parse_gate_application()
        elif token.type == 'ID' and peek.type == 'EQUALS':
             # Assignment: bit = measure q;
             target_bit = Identifier(self.eat('ID').value)
             self.eat('EQUALS')
             if self.current_token and self.current_token.type == 'MEASURE':
                 measure_stmt = self.parse_measure_statement()
                 if measure_stmt:
                    measure_stmt.classical_bit = target_bit # Override
                 return measure_stmt
             else:
                 current_type = self.current_token.type if self.current_token else 'EOF'
                 current_line = self.current_token.line if self.current_token else 0
                 current_col = self.current_token.column if self.current_token else 0
                 raise CompilationError(f"Expected MEASURE after assignment but got {current_type}", current_line, current_col)
        elif token.type == 'MEASURE':
            return self.parse_measure_statement()
        elif token.type == 'APPLY':
             stmt = self.parse_stimulus_application()
             self.eat('SEMICOLON')
             return stmt
        elif token.type in ['H', 'X', 'Y', 'Z', 'CNOT', 'CX']:
             return self.parse_quantum_gate()
        else:
            raise CompilationError(f"Invalid statement in circuit: {token.type}", token.line, token.column)

    def parse_qubit_declaration(self):
        name_token = self.eat('ID')
        name = Identifier(name_token.value)
        self.eat('COLON')
        type_token = self.eat('ID')
        self.eat('SEMICOLON')
        
        if type_token.value == 'qubit':
            return QubitDeclaration(name=name, line=name_token.line, column=name_token.column)
        elif type_token.value == 'bit':
            return ClassicalBitDefinition(name=name, line=name_token.line, column=name_token.column)
        else:
            raise CompilationError(f"Expected 'qubit' or 'bit' but got '{type_token.value}'", type_token.line, type_token.column)

    def parse_gate_application(self):
        gate_token = self.current_token
        if not gate_token:
            raise CompilationError("Unexpected end of file while parsing gate application.", 1, 1)
        gate_name = Identifier(self.eat('ID').value)
        self.eat('LPAREN')
        arguments = self.parse_arguments()
        self.eat('RPAREN')
        self.eat('ON')
        qubits = self.parse_qubit_list()
        self.eat('SEMICOLON')
        return GateApplication(gate_name=gate_name, arguments=arguments, qubits=qubits, line=gate_token.line, column=gate_token.column)

    def parse_measure_statement(self):
        measure_token = self.current_token
        if not measure_token:
            raise CompilationError("Unexpected end of file while parsing measure statement.", 1, 1)
        self.eat('MEASURE')
        qubit = Identifier(self.eat('ID').value)
        
        target_bit = None
        # Handle ON syntax: MEASURE q1 ON c1;
        if self.current_token and self.current_token.type == 'ON':
            self.eat('ON')
            target_bit = Identifier(self.eat('ID').value)
        
        # Handle arrow syntax: measure q -> m;
        elif self.current_token and self.current_token.type == 'ARROW':
            self.eat('ARROW')
            target_bit = Identifier(self.eat('ID').value)
        
        # Handle traditional syntax: measure q to m;
        elif self.current_token and self.current_token.type == 'TO':
            self.eat('TO')
            target_bit = Identifier(self.eat('ID').value)
            
        elif self.current_token and self.current_token.type == 'SEMICOLON':
            # Allow measure q; (assignment handled elsewhere or implicit)
            pass
        
        else:
            self.error_reporter.report("Expected 'ON', '->', 'TO' or ';' after qubit in measure statement", measure_token.line, measure_token.column)
            return None
            
        self.eat('SEMICOLON')
        return MeasureStatement(qubit=qubit, classical_bit=target_bit, line=measure_token.line, column=measure_token.column)

    def parse_arguments(self):
        args = []
        if self.current_token and self.current_token.type != 'RPAREN':
            while self.current_token and self.current_token.type != 'RPAREN':
                args.append(Argument(value=self.eat('ID').value))
                if self.current_token and self.current_token.type == 'COMMA':
                    self.eat('COMMA')
                else:
                    break
        return args

    def parse_quantum_operation(self):
        token = self.current_token
        if not token:
            raise CompilationError("Unexpected end of file in quantum operation", 1, 1)
        if token.type in ['H', 'X', 'Y', 'Z', 'CNOT']:
            name = Identifier(self.eat(token.type).value)
        elif token.type == 'ID':
            name = Identifier(self.eat('ID').value)
        else:
            raise CompilationError(f"Unexpected token {token.type} in quantum operation", token.line, token.column)

        self.eat('LPAREN')
        qubit_refs = self.parse_identifier_list()
        self.eat('RPAREN')
        self.eat('SEMICOLON')
        return QuantumOperation(name=name, qubit_refs=qubit_refs)

    
