import re
import json
import sys

def parse_genbank(file_path):
    features = []
    locus_size = 0
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    in_features = False
    current_feature = None
    
    for line in lines:
        line = line.strip()
        
        if line.startswith("LOCUS"):
            parts = line.split()
            # Usually size is the 3rd element, e.g. "LOCUS HAWRA 17243 bp"
            for part in parts:
                if part.isdigit():
                    locus_size = int(part)
                    break
        
        if line.startswith("FEATURES"):
            in_features = True
            continue
            
        if line.startswith("ORIGIN"):
            in_features = False
            continue
            
        if in_features:
            # Check for feature start (e.g., "     gene            801..2800")
            # Indentation matters in GenBank but we stripped it. 
            # Real GenBank feature keys start at col 5.
            # Let's use regex for safety.
            
            # Match feature type and location
            match = re.match(r'^([a-zA-Z_]+)\s+(\d+)\.\.(\d+)', line)
            if match:
                if current_feature:
                    features.append(current_feature)
                
                ftype = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                current_feature = {
                    "type": ftype,
                    "start": start,
                    "end": end,
                    "label": ftype,
                    "note": ""
                }
            elif current_feature:
                # Parse qualifiers like /gene="psaA" or /note="..."
                if line.startswith('/'):
                    q_match = re.match(r'/([^=]+)="?([^"]*)"?', line)
                    if q_match:
                        key = q_match.group(1)
                        val = q_match.group(2)
                        if key == "gene" or key == "label":
                            current_feature["label"] = val
                        elif key == "note":
                            current_feature["note"] += val + " "

    if current_feature:
        features.append(current_feature)
        
    return {"size": locus_size, "features": features}

files = [
    "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/plasmids/HAWRA_FINAL_VALIDATED.gb",
    "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/plasmids/HAWRA_PLASMID_v3.gb"
]

results = {}
for fp in files:
    name = fp.split('/')[-1]
    try:
        results[name] = parse_genbank(fp)
    except Exception as e:
        results[name] = {"error": str(e)}

print(json.dumps(results, indent=2))
