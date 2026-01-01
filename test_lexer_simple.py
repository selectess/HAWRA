from arbol.compiler.lexer import Lexer
import sys

source = """// HAWRA "First Bloom" Experiment
genes {
    gene "psaA"
}"""

lexer = Lexer(source)
for token in lexer.tokens():
    print(token)
