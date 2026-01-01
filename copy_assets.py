import shutil
import os

# Source directories to check for assets
source_dirs = [
    "/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE/01_manuscript/figures",
    "/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE/04_multimedia",
    "/Users/mehdiwhb/Desktop/HAWRA/HAWRA_SUBMISSION_BUNDLE_2025_12_21/01_Manuscript/figures",
    "/Users/mehdiwhb/Desktop/HAWRA/06_publication/latex/figures",
    "/Users/mehdiwhb/Desktop/HAWRA/05_data/results"
]

target_dir = "/Users/mehdiwhb/Desktop/HAWRA/HAWRA_FINAL_SUBMISSION_2025_12_21/01_Manuscript/figures"
multimedia_target = "/Users/mehdiwhb/Desktop/HAWRA/HAWRA_FINAL_SUBMISSION_2025_12_21/04_Multimedia"

os.makedirs(target_dir, exist_ok=True)
os.makedirs(multimedia_target, exist_ok=True)

extensions = ('.png', '.jpg', '.jpeg', '.gif', '.pdf', '.svg')

copied_count = 0

for s_dir in source_dirs:
    if os.path.exists(s_dir):
        print(f"Checking {s_dir}...")
        for item in os.listdir(s_dir):
            if item.lower().endswith(extensions):
                s = os.path.join(s_dir, item)
                d = os.path.join(target_dir, item)
                if not os.path.exists(d):
                    try:
                        shutil.copy2(s, d)
                        copied_count += 1
                        # Also copy to multimedia if it's from a multimedia source or is a gif
                        if "multimedia" in s_dir.lower() or item.lower().endswith('.gif'):
                            shutil.copy2(s, os.path.join(multimedia_target, item))
                    except Exception as e:
                        print(f"Error copying {item}: {e}")

print(f"Copy completed. {copied_count} files copied.")
