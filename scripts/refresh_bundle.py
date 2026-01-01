import os
import shutil
from datetime import datetime

# Configuration
PROJECT_ROOT = "/Users/mehdiwhb/Desktop/HAWRA"
BUNDLE_NAME = f"HAWRA_SUBMISSION_BUNDLE_{datetime.now().strftime('%Y_%m_%d')}"
BUNDLE_PATH = os.path.join(PROJECT_ROOT, BUNDLE_NAME)

STRUCTURE = {
    "00_Manifesto": [
        "06_publication/HAWRA_Manifesto.md",
        "06_publication/HAWRA_Nature_Submission.md"
    ],
    "01_Manuscript": [
        "06_publication/latex/main.tex",
        "06_publication/latex/references.bib",
        "06_publication/latex/sections",
        "06_publication/latex/figures"
    ],
    "02_Code": [
        "arbol",
            "bioos",
            "02_arbol_interface/jetson_client"
    ],
    "03_Data": [
        "06_publication/latex/figures/phytoqmmml_convergence.json"
    ],
    "04_Multimedia": [
        "06_publication/latex/figures"
    ],
    "05_Supplementary": [
        "00_docs/scientific/SOP_Implementation_Guide.md",
        "00_docs/project_management/HAWRA_Scientific_FAQ_Defense.md"
    ]
}

def create_bundle():
    if os.path.exists(BUNDLE_PATH):
        shutil.rmtree(BUNDLE_PATH)
    
    os.makedirs(BUNDLE_PATH)
    
    for folder, items in STRUCTURE.items():
        folder_path = os.path.join(BUNDLE_PATH, folder)
        os.makedirs(folder_path, exist_ok=True)
        
        for item in items:
            src = os.path.join(PROJECT_ROOT, item)
            if not os.path.exists(src):
                print(f"Warning: {src} does not exist.")
                continue
            
            dest = os.path.join(folder_path, os.path.basename(item))
            if os.path.isdir(src):
                shutil.copytree(src, dest, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dest)
    
    # Create README in bundle
    readme_content = f"""# HAWRA Submission Bundle
Date: {datetime.now().strftime('%Y-%m-%d')}
Conceived & Engineered by: Mehdi Wahbi (Director, Move37 Initiative)
Implemented by: Move37 AI Team (Move37 Initiative)
Contact: m.wahbi.move37@atomicmail.io
DOI: 10.5281/zenodo.17908061

## Contents
- 00_Manifesto: Core vision and submission letters.
- 01_Manuscript: Full LaTeX source for the scientific paper.
- 02_Code: Implementation of ARBOL, BioOS, and Jetson Interface.
- 03_Data: Simulation results and validation data.
- 04_Multimedia: High-resolution figures and diagrams.
- 05_Supplementary: SOPs and technical FAQs.
"""
    with open(os.path.join(BUNDLE_PATH, "README.md"), "w") as f:
        f.write(readme_content)
    
    print(f"Bundle created at: {BUNDLE_PATH}")

if __name__ == "__main__":
    create_bundle()
