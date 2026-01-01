import shutil
import os

src_dir = '/Users/mehdiwhb/Desktop/HAWRA/06_publication/latex/figures/'
dst_dir = '/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE_2025_12_21/01_Manuscript/figures/'

if not os.path.exists(dst_dir):
    os.makedirs(dst_dir)

files = os.listdir(src_dir)
for f in files:
    src_path = os.path.join(src_dir, f)
    dst_path = os.path.join(dst_dir, f)
    if os.path.isfile(src_path):
        try:
            shutil.copy2(src_path, dst_path)
            print(f"Copied {f}")
        except Exception as e:
            print(f"Failed to copy {f}: {e}")
