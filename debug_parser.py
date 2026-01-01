#!/usr/bin/env python3

from arbol.compiler.lexer import Lexer
from arbol.compiler.parser import Parser
from arbol.compiler.error import ErrorReporter

def debug_parse():
    with open('arbol/test.arbol', 'r') as f:
        arbol_code = f.read()

    error_reporter = ErrorReporter()
    lexer = Lexer(arbol_code)
    tokens = list(lexer.tokens())
    
    print("Tokens:")
    for i, token in enumerate(tokens):
        print(f"{i:2d}: Line {token.line:2d}, Col {token.column:2d} - {token.type:15s} '{token.value}'")
    
    print("\nParsing...")
    parser = Parser(lexer, error_reporter)
    
    # Manually step through parsing to see what happens
    pos = 0
    statements = []
    
    while pos < len(tokens):
        token = tokens[pos]
        print(f"\nProcessing token at pos {pos}: {token}")
        
        if token.type == 'STIMULUS':
            print("  -> Parsing stimulus definition")
        elif token.type == 'LOGICAL_QUBIT':
            print("  -> Parsing logical qubit definition")
        elif token.type == 'CLASSICAL_BIT':
            print("  -> Parsing classical bit definition")
        elif token.type == 'GATE':
            print("  -> Parsing gate definition")
        elif token.type == 'H':
            print("  -> Parsing H gate application")
        elif token.type == 'MEASURE':
            print("  -> Parsing measure statement")
        else:
            print(f"  -> Unknown token type: {token.type}")
            
        pos += 1
        if pos > 35:  # Limit for debugging
            break
    
    print(f"\nErrors: {len(error_reporter.get_errors())}")
    for error in error_reporter.get_errors():
        print(f"  - {error}")

if __name__ == "__main__":
    debug_parse()