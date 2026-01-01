import os
import json
import hashlib

def calculate_sha256(file_path):
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

bundle_dir = "HAWRA_SUBMISSION_BUNDLE_2025_12_21"
manifest = {
    "project": "HAWRA",
    "version": "v1.0-submission",
    "date": "2025-12-21",
    "author": "Mehdi Wahbi",
    "files": []
}

if os.path.exists(bundle_dir):
    for root, dirs, files in os.walk(bundle_dir):
        for file in files:
            file_path = os.path.join(root, file)
            rel_path = os.path.relpath(file_path, bundle_dir)
            manifest["files"].append({
                "path": rel_path,
                "sha256": calculate_sha256(file_path)
            })

    with open(os.path.join(bundle_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=4)
    print("Manifest generated successfully.")
else:
    print("Bundle directory not found.")
