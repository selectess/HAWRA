import os
import re

def escape_underscores(content):
    # Split by lstlisting blocks
    parts = re.split(r'(\\begin\{lstlisting\}.*?\\end\{lstlisting\})', content, flags=re.DOTALL)
    for i in range(len(parts)):
        if not parts[i].startswith('\\begin{lstlisting}'):
            # Escape underscores not preceded by backslash
            parts[i] = re.sub(r'(?<!\\)_', r'\\_', parts[i])
    return "".join(parts)

directory = "/Users/mehdiwhb/Desktop/HAWRA/HAWRA_DOCUMENTARY_REPORT/sections"
for filename in os.listdir(directory):
    if filename.endswith(".tex"):
        path = os.path.join(directory, filename)
        with open(path, 'r') as f:
            content = f.read()
        new_content = escape_underscores(content)
        with open(path, 'w') as f:
            f.write(new_content)
