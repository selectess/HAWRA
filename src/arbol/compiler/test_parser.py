import unittest
from compiler.lexer import Lexer
from compiler.parser import Parser
from compiler.arbol_ast import (
    Program, CircuitDefinition, GateDefinition, QubitDeclaration, 
    GateApplication, MeasureStatement, Identifier, Argument
)
from compiler.error import ErrorReporter

class TestParser(unittest.TestCase):
    def _parse_code(self, code):
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()
        return ast, error_reporter

    def test_circuit_definition_with_qubit_declaration(self):
        code = """
        CIRCUIT test() {
            q1: qubit;
        }
        """
        ast, error_reporter = self._parse_code(code)
        self.assertFalse(error_reporter.has_errors())
        expected_ast = Program(statements=[
            CircuitDefinition(
                name=Identifier('test'), 
                parameters=[], 
                body=[
                    QubitDeclaration(name=Identifier('q1'))
                ]
            )
        ])
        self.assertEqual(ast, expected_ast)

    def test_gate_definition(self):
        code = """
        GATE H(q) {
            X(q);
        }
        """
        ast, error_reporter = self._parse_code(code)
        self.assertFalse(error_reporter.has_errors())
        # Note: La validation complète de la sémantique du corps de la porte est du ressort du compilateur.
        # Le parser vérifie simplement la structure.
        self.assertIsInstance(ast.statements[0], GateDefinition)
        self.assertEqual(ast.statements[0].name, Identifier('H'))

if __name__ == '__main__':
    unittest.main()
