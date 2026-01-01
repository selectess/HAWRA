import os
import shutil

bundle_root = 'HAWRA_SUBMISSION_BUNDLE_2025_12_21'
dirs = [
    '00_Manifesto',
    '01_Manuscript',
    '01_Manuscript/sections',
    '02_Code',
    '02_Code/bioos',
    '02_Code/simulator',
    '02_Code/arbol',
    '03_Data',
    '04_Multimedia',
    '05_Supplementary'
]

# Create directories
for d in dirs:
    path = os.path.join(bundle_root, d)
    os.makedirs(path, exist_ok=True)
    print(f"Created {path}")

# Helper to copy files
def copy_safe(src, dst):
    try:
        if os.path.isdir(src):
            shutil.copytree(src, dst, dirs_exist_ok=True)
        else:
            shutil.copy2(src, dst)
        print(f"Copied {src} to {dst}")
    except Exception as e:
        print(f"Error copying {src}: {e}")

# 00_Manifesto
copy_safe('06_publication/HAWRA_Manifesto.md', os.path.join(bundle_root, '00_Manifesto/'))

# 01_Manuscript
copy_safe('06_publication/latex/main.pdf', os.path.join(bundle_root, '01_Manuscript/'))
copy_safe('06_publication/latex/main.tex', os.path.join(bundle_root, '01_Manuscript/'))
copy_safe('06_publication/latex/references.bib', os.path.join(bundle_root, '01_Manuscript/'))
copy_safe('06_publication/latex/sections/', os.path.join(bundle_root, '01_Manuscript/sections/'))
copy_safe('06_publication/HAWRA_Nature_Submission.md', os.path.join(bundle_root, '01_Manuscript/'))
copy_safe('06_publication/HAWRA_Full_Scientific_Paper.md', os.path.join(bundle_root, '01_Manuscript/'))

# 02_Code
copy_safe('bioos/', os.path.join(bundle_root, '02_Code/bioos/'))
copy_safe('03_unified_simulator/src/', os.path.join(bundle_root, '02_Code/simulator/'))
copy_safe('02_arbol_interface/', os.path.join(bundle_root, '02_Code/arbol/'))

# 03_Data
copy_safe('01_genomics/plasmids/validated/HAWRA_FINAL_VALIDATED.gb', os.path.join(bundle_root, '03_Data/'))
copy_safe('03_unified_simulator/results/', os.path.join(bundle_root, '03_Data/'))

# 04_Multimedia
copy_safe('06_publication/latex/figures/', os.path.join(bundle_root, '04_Multimedia/'))

# 05_Supplementary
copy_safe('00_docs/concepts/SOP_ficus_regeneration.md', os.path.join(bundle_root, '05_Supplementary/'))
copy_safe('00_docs/concepts/PQPE_protocol.md', os.path.join(bundle_root, '05_Supplementary/'))
copy_safe('00_docs/analysis/sensitivity_analysis_report.md', os.path.join(bundle_root, '05_Supplementary/'))

print("Bundle finalization complete.")
