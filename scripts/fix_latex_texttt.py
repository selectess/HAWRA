import os
import re

directory = '/Users/mehdiwhb/Desktop/HAWRA/06_publication/latex/sections/'

def fix_texttt(match):
    full_match = match.group(0)
    # Only replace if there's an underscore
    if '_' in full_match:
        fixed = full_match.replace('_', '\\_')
        return fixed
    return full_match

# Pattern to match \texttt{...}
# We use a non-greedy match for the content
texttt_pattern = re.compile(r'\\texttt\{.*?\}')

files_fixed = 0
for filename in os.listdir(directory):
    if filename.endswith('.tex'):
        filepath = os.path.join(directory, filename)
        with open(filepath, 'r') as f:
            content = f.read()
        
        new_content = texttt_pattern.sub(fix_texttt, content)
        
        if new_content != content:
            with open(filepath, 'w') as f:
                f.write(new_content)
            print(f'Fixed {filename}')
            files_fixed += 1

print(f'Total files fixed: {files_fixed}')
