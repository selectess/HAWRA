import sys
sys.path.append('/Users/mehdiwhb/Desktop/HAWRA')

def test_simple():
    assert 1 == 1

def test_import():
    from arbol.compiler.lexer import Lexer
    assert Lexer is not None