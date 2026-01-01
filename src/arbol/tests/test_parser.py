import unittest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from compiler.lexer import Lexer
from compiler.parser import Parser
from compiler.arbol_ast import *

from compiler.error import ErrorReporter

class TestParser(unittest.TestCase):

    def test_empty_program(self):
        code = ""
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()
        self.assertEqual(ast, Program(statements=[]))

    def test_circuit_definition(self):
        code = "CIRCUIT my_circuit { }"
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()
        expected_ast = Program(statements=[
            CircuitDefinition(name=Identifier(name='my_circuit'), parameters=[], body=[])
        ])
        self.assertEqual(ast, expected_ast)

    def test_logical_qubit_definition(self):
        code = "LOGICAL_QUBIT q1 [2] IS { activator: p1, repressor: p2 };"
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()

        self.assertTrue(ast.statements, "Le parser n'a produit aucun statement.")

        actual_line = ast.statements[0].line
        actual_column = ast.statements[0].column

        expected_ast = Program(statements=[
            LogicalQubitDefinition(
                name=Identifier(name='q1'),
                size=2,
                is_clause=IsClause(promoters=[
                    Promoter(type='activator', name='p1'),
                    Promoter(type='repressor', name='p2')
                ]),
                line=actual_line,
                column=actual_column
            )
        ])
        self.assertEqual(ast, expected_ast)

    def test_measure_statement_on(self):
        code = "MEASURE q1 ON c1;"
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()
        if ast.statements:
            actual_line = ast.statements[0].line
            actual_column = ast.statements[0].column
        expected = Program(statements=[
            MeasureStatement(qubit=Identifier(name='q1'), classical_bit=Identifier(name='c1'), line=actual_line, column=actual_column)
        ])
        self.assertEqual(ast, expected)

    def test_measure_statement_arrow(self):
        code = "MEASURE q1 -> m1;"
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()
        if ast.statements:
            actual_line = ast.statements[0].line
            actual_column = ast.statements[0].column
        expected = Program(statements=[
            MeasureStatement(qubit=Identifier(name='q1'), classical_bit=Identifier(name='m1'), line=actual_line, column=actual_column)
        ])
        self.assertEqual(ast, expected)

    def test_measure_statement_to(self):
        code = "MEASURE q1 TO m1;"
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()
        if ast.statements:
            actual_line = ast.statements[0].line
            actual_column = ast.statements[0].column
        expected = Program(statements=[
            MeasureStatement(qubit=Identifier(name='q1'), classical_bit=Identifier(name='m1'), line=actual_line, column=actual_column)
        ])
        self.assertEqual(ast, expected)

    def test_single_gate_on_qubit(self):
        code = "H on q1;"
        lexer = Lexer(code)
        error_reporter = ErrorReporter()
        parser = Parser(lexer, error_reporter)
        ast = parser.parse()
        stmt = ast.statements[0]
        # Vérifie type et contenu minimal
        self.assertIsInstance(stmt, GateApplication)
        self.assertEqual(stmt.gate_name, Identifier(name='H'))
        self.assertEqual(stmt.qubits, [Identifier(name='q1')])

if __name__ == '__main__':
    unittest.main()
