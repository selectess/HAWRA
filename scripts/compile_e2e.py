import json
import sys
from arbol.compiler.lexer import Lexer
from arbol.compiler.parser import Parser
from arbol.compiler.compiler import Compiler
from arbol.compiler.error import ErrorReporter

def main():
    """Compile un fichier Arbol en .bsim.json."""
    if len(sys.argv) != 2:
        print("Usage: python compile_e2e.py <fichier_arbol>")
        sys.exit(1)

    print("Lancement de la compilation...")
    filepath = sys.argv[1]
    with open(filepath, 'r') as f:
        source_code = f.read()
    print("Source code lu.")

    error_reporter = ErrorReporter()
    lexer = Lexer(source_code)
    print("Lexer créé.")

    parser = Parser(lexer, error_reporter)
    print("Parser créé.")
    ast = parser.parse()
    print("AST généré.")

    if error_reporter.has_errors():
        error_reporter.display()
        sys.exit(1)

    compiler = Compiler(error_reporter)
    print("Compilateur créé.")
    assembly = compiler.compile(ast)
    print("Assembly généré.")

    if error_reporter.has_errors():
        error_reporter.display()
        sys.exit(1)

    output_path = filepath.replace('.arbol', '.bsim.json')
    with open(output_path, 'w') as f:
        json.dump(assembly, f, indent=4)

    print(f"Compilation réussie : {output_path}")

if __name__ == "__main__":
    main()
