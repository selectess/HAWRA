import shutil
import os

def copy_tree(src, dst):
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            if item != '__pycache__':
                copy_tree(s, d)
        else:
            shutil.copy2(s, d)
            print(f"Copied {item}")

# 1. Copier le code
src_code = '/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE/02_codebase/'
dst_code = '/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE_2025_12_21/02_Code/'
copy_tree(src_code, dst_code)

# 2. Copier les figures
src_figs = '/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE/01_manuscript/figures/'
dst_figs = '/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE_2025_12_21/01_Manuscript/figures/'
copy_tree(src_figs, dst_figs)

print("Migration complete.")
