import sys
import pkg_resources

REQUIRED_PACKAGES = [
    'numpy',
    'biopython',
    'matplotlib',
    'pandas',
    'plotly',
    'qutip',
    'kaleido',
    'jupyter'
]

print(f"Python Version: {sys.version}")
print("-" * 20)

for package in REQUIRED_PACKAGES:
    try:
        dist = pkg_resources.get_distribution(package)
        print(f"[OK] {dist.key} ({dist.version})")
    except pkg_resources.DistributionNotFound:
        print(f"[FAIL] {package} is not installed.")
    except Exception as e:
        print(f"[ERROR] An error occurred while checking {package}: {e}")

print("-" * 20)
print("Environment check complete.")
