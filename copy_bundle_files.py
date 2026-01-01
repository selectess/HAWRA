import os
import shutil

src_sections = '/Users/mehdiwhb/Desktop/HAWRA/06_publication/latex/sections/'
dst_sections = '/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE_2025_12_21/01_Manuscript/sections/'

src_figures = '/Users/mehdiwhb/Desktop/HAWRA/06_publication/latex/figures/'
dst_figures = '/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE_2025_12_21/01_Manuscript/figures/'

def copy_dir(src, dst):
    if not os.path.exists(dst):
        os.makedirs(dst)
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isfile(s):
            shutil.copy2(s, d)
            print(f"Copied {s} to {d}")

try:
    copy_dir(src_sections, dst_sections)
    copy_dir(src_figures, dst_figures)
    print("Success")
except Exception as e:
    print(f"Error: {e}")
