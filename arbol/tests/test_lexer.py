import unittest
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from compiler.lexer import Lexer, Token

class TestLexer(unittest.TestCase):

    def test_keywords(self):
        code = "CIRCUIT APPLY LOGICAL_QUBIT CLASSICAL_BIT GATE STIMULUS IS AT SITE ON MEASURE H X Y Z CNOT ACTIVATOR REPRESSOR TO WITH"
        lexer = Lexer(code)
        tokens = list(lexer.tokens())
        expected_tokens = [
            Token('CIRCUIT', 'CIRCUIT', 1, 1),
            Token('APPLY', 'APPLY', 1, 9),
            Token('LOGICAL_QUBIT', 'LOGICAL_QUBIT', 1, 15),
            Token('CLASSICAL_BIT', 'CLASSICAL_BIT', 1, 29),
            Token('GATE', 'GATE', 1, 43),
            Token('STIMULUS', 'STIMULUS', 1, 48),
            Token('IS', 'IS', 1, 57),
            Token('AT', 'AT', 1, 60),
            Token('SITE', 'SITE', 1, 63),
            Token('ON', 'ON', 1, 68),
            Token('MEASURE', 'MEASURE', 1, 71),
            Token('H', 'H', 1, 79),
            Token('X', 'X', 1, 81),
            Token('Y', 'Y', 1, 83),
            Token('Z', 'Z', 1, 85),
            Token('CNOT', 'CNOT', 1, 87),
            Token('ACTIVATOR', 'ACTIVATOR', 1, 92),
            Token('REPRESSOR', 'REPRESSOR', 1, 102),
            Token('TO', 'TO', 1, 112),
            Token('WITH', 'WITH', 1, 115),
        ]
        self.assertEqual(tokens, expected_tokens)

    def test_identifiers(self):
        code = "my_circuit my_qubit my_gate"
        lexer = Lexer(code)
        tokens = list(lexer.tokens())
        expected_tokens = [
            Token('ID', 'my_circuit', 1, 1),
            Token('ID', 'my_qubit', 1, 12),
            Token('ID', 'my_gate', 1, 21),
        ]
        self.assertEqual(tokens, expected_tokens)

    def test_numbers(self):
        code = "123 45.67 0.89"
        lexer = Lexer(code)
        tokens = list(lexer.tokens())
        expected_tokens = [
            Token('NUMBER', '123', 1, 1),
            Token('NUMBER', '45.67', 1, 5),
            Token('NUMBER', '0.89', 1, 11),
        ]
        self.assertEqual(tokens, expected_tokens)

    def test_strings(self):
        code = '"hello world" "another string"'
        lexer = Lexer(code)
        tokens = list(lexer.tokens())
        expected_tokens = [
            Token('STRING', 'hello world', 1, 1),
            Token('STRING', 'another string', 1, 15),
        ]
        self.assertEqual(tokens, expected_tokens)

    def test_operators_and_delimiters(self):
        code = "{ } ( ) [ ] , : ; ->"
        lexer = Lexer(code)
        tokens = list(lexer.tokens())
        expected_tokens = [
            Token('LBRACE', '{', 1, 1),
            Token('RBRACE', '}', 1, 3),
            Token('LPAREN', '(', 1, 5),
            Token('RPAREN', ')', 1, 7),
            Token('LBRACKET', '[', 1, 9),
            Token('RBRACKET', ']', 1, 11),
            Token('COMMA', ',', 1, 13),
            Token('COLON', ':', 1, 15),
            Token('SEMICOLON', ';', 1, 17),
            Token('ARROW', '->', 1, 19),
        ]
        self.assertEqual(tokens, expected_tokens)

    def test_comments_and_whitespace(self):
        code = """
        (* this is a comment *)
        CIRCUIT my_circuit {
            (* another comment *)
        }
        """
        lexer = Lexer(code)
        tokens = list(lexer.tokens())
        expected_tokens = [
            Token('CIRCUIT', 'CIRCUIT', 3, 9),
            Token('ID', 'my_circuit', 3, 17),
            Token('LBRACE', '{', 3, 28),
            Token('RBRACE', '}', 5, 9),
        ]
        self.assertEqual(tokens, expected_tokens)

    def test_unknown_token(self):
        code = "@"
        lexer = Lexer(code)
        with self.assertRaises(Exception):
            list(lexer.tokens())

if __name__ == '__main__':
    unittest.main()
