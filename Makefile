PY=python3

pdf:
	$(PY) scripts/build_manuscript_pdf.py

preprint:
	$(PY) scripts/build_preprint_bundle.py

full:
	$(PY) scripts/build_full_archive.py

test:
	$(PY) -m unittest discover -s 03_unified_simulator/tests -p "test_*.py"

ci: pdf test preprint