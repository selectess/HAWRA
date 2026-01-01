import re

class Token:
    def __init__(self, type, value, line, column):
        self.type = type
        self.value = value
        self.line = line
        self.column = column

    def __eq__(self, other):
        return (self.type == other.type and
                self.value == other.value and
                self.line == other.line and
                self.column == other.column)

    def __repr__(self):
        return f"Token({self.type}, {self.value!r}, {self.line}, {self.column})"

class Lexer:
    def __init__(self, text):
        self.text = text
        self.pos = 0
        self.line = 1
        self.column = 1
        self.token_specification = [
            ('COMMENT_SINGLE', r'//[^\n]*'),       # Single line comments
            ('COMMENT',     r'\(\*.*?\*\)'),  # Comments (* ... *)
            ('STRING',      r'"[^"]*"'),           # String literals
            ('NUMBER',      r'\d+(\.\d*)?([eE][+-]?\d+)?'), # Integer, decimal, or scientific notation
            ('CIRCUIT',     r'circuit\b'),             # Keywords
            ('APPLY',       r'apply\b'),
            ('TO',          r'to\b'),
            ('WITH',        r'with\b'),
            ('ACTIVATOR',   r'activator\b'),
            ('REPRESSOR',   r'repressor\b'),
            
            ('STIMULUS',    r'stimulus\b'),
            ('LOGICAL_QUBIT', r'LOGICAL_QUBIT\b'),
            ('CLASSICAL_BIT', r'classical_bit\b'),
            ('IS',          r'is\b'),
            ('AT',          r'at\b'),
            ('SITE',        r'site\b'),
            ('GATE',        r'gate\b'),
            ('ON',          r'on\b'),
            ('OVER',        r'over\b'),
            ('CONFIG',      r'config\b'),
            ('RUN',         r'run\b'),
            ('GENES',       r'genes\b'),
            ('GRN',         r'grn\b'),
            ('ACTIVATES',   r'activates\b'),
            ('REPRESSES',   r'represses\b'),
            ('BASAL_RATE',  r'basal_rate\b'),
            ('DEGRADATION_RATE', r'degradation_rate\b'),
            ('WEIGHT',      r'weight\b'),
            ('HILL_COEFFICIENT', r'hill_coefficient\b'),
            ('HALF_MAX_CONCENTRATION', r'half_max_concentration\b'),
            ('STEP',        r'step\b'),
            ('H',           r'H\b'),
            ('X',           r'X\b'),
            ('Y',           r'Y\b'),
            ('Z',           r'Z\b'),
            ('CNOT',        r'CNOT\b'),
            ('MEASURE',     r'measure\b'),
            ('ARROW',       r'->'),                     # Operators and delimiters
            ('LBRACE',      r'\{'),
            ('RBRACE',      r'\}'),
            ('LPAREN',      r'\('),
            ('RPAREN',      r'\)'),
            ('LBRACKET',    r'\['),
            ('RBRACKET',    r'\]'),
            ('COLON',       r':'),
            ('EQUALS',      r'='),
            ('COMMA',       r','),
            ('SEMICOLON',   r';'),
            ('ID',          r'[A-Za-z_][A-Za-z0-9_]*'), # Identifiers
            ('NEWLINE',     r'\n'),                     # Line endings
            ('SKIP',        r'[ \t]+'),                 # Skip whitespace
            ('MISMATCH',    r'.'),                      # Any other character
        ]
        self.tok_regex = re.compile('|'.join('(?P<%s>%s)' % pair for pair in self.token_specification), re.IGNORECASE | re.DOTALL)

    def tokens(self):
        for mo in self.tok_regex.finditer(self.text):
            kind = mo.lastgroup
            value = mo.group()
            column = self.column
            if kind == 'NEWLINE':
                self.line += 1
                self.column = 1
                continue
            elif kind == 'SKIP':
                self.column += len(value)
                continue
            elif kind == 'COMMENT' or kind == 'COMMENT_SINGLE':
                newlines = value.count('\n')
                if newlines > 0:
                    self.line += newlines
                    self.column = len(value) - value.rfind('\n')
                else:
                    self.column += len(value)
                continue
            elif kind == 'STRING':
                original_len = len(value)
                value = value[1:-1] # Remove quotes
                token = Token(kind, value, self.line, column)
                self.column += original_len
                yield token
                continue
            elif kind == 'MISMATCH':
                raise RuntimeError(f'{value!r} unexpected on line {self.line}')
            
            if kind == 'MISMATCH':
                raise RuntimeError(f'{value!r} unexpected on line {self.line}')
            
            token = Token(kind, value, self.line, column)
            self.column += len(value)
            yield token

if __name__ == '__main__':
    code = '''
    (* Exemple de code ARBOL : Contrôle lumineux d'un gène *)

    circuit light_switch {
        CRY2 -> P700 [ activator: light_blue ]
    }

    apply light to plant with {
        wavelength = "488nm",
        intensity = "high",
        duration = "60min"
    }
    '''
    
    lexer = Lexer(code)
    for token in lexer.tokens():
        print(token)
