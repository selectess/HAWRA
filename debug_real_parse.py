#!/usr/bin/env python3

from arbol.compiler.lexer import Lexer
from arbol.compiler.parser import Parser
from arbol.compiler.error import ErrorReporter

def debug_real_parse():
    with open('arbol/test.arbol', 'r') as f:
        arbol_code = f.read()

    error_reporter = ErrorReporter()
    lexer = Lexer(arbol_code)
    parser = Parser(lexer, error_reporter)
    
    # Parse and see what happens
    try:
        ast = parser.parse()
        print("Parsing completed successfully!")
        print(f"AST statements: {len(ast.statements)}")
        for i, stmt in enumerate(ast.statements):
            print(f"  {i}: {type(stmt).__name__}")
            if hasattr(stmt, 'name'):
                print(f"      name: {stmt.name}")
            if hasattr(stmt, 'line'):
                print(f"      line: {stmt.line}")
    except Exception as e:
        print(f"Exception during parsing: {e}")
        import traceback
        traceback.print_exc()
    
    print(f"\nErrors: {len(error_reporter.get_errors())}")
    for error in error_reporter.get_errors():
        print(f"  - {error}")

if __name__ == "__main__":
    debug_real_parse()