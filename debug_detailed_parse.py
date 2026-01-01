#!/usr/bin/env python3

from arbol.compiler.lexer import Lexer
from arbol.compiler.parser import Parser
from arbol.compiler.error import ErrorReporter
import traceback

def debug_detailed_parse():
    with open('arbol/test.arbol', 'r') as f:
        arbol_code = f.read()

    error_reporter = ErrorReporter()
    lexer = Lexer(arbol_code)
    
    # Monkey patch the parse_statement method to add debug output
    original_parse_statement = Parser.parse_statement
    
    def debug_parse_statement(self):
        token = self.current_token
        if not token:
            return None
        
        print(f"parse_statement called with token: {token}")
        
        try:
            result = original_parse_statement(self)
            if result:
                print(f"  -> Successfully parsed: {type(result).__name__}")
            else:
                print(f"  -> No result")
            return result
        except Exception as e:
            print(f"  -> Exception: {e}")
            traceback.print_exc()
            raise
    
    Parser.parse_statement = debug_parse_statement
    
    parser = Parser(lexer, error_reporter)
    
    try:
        ast = parser.parse()
        print("Parsing completed successfully!")
        print(f"AST statements: {len(ast.statements)}")
        for i, stmt in enumerate(ast.statements):
            print(f"  {i}: {type(stmt).__name__}")
            if hasattr(stmt, 'line') and stmt.line is not None:
                print(f"      line: {stmt.line}")
    except Exception as e:
        print(f"Exception during parsing: {e}")
        traceback.print_exc()
    
    print(f"\nErrors: {len(error_reporter.get_errors())}")
    for error in error_reporter.get_errors():
        print(f"  - {error}")

if __name__ == "__main__":
    debug_detailed_parse()