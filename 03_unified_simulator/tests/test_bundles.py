import os
import unittest
import zipfile

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
PREPRINT_ZIP = os.path.join(ROOT, 'dist', 'preprint_submission', 'HAWRA_preprint_bundle.zip')
FULL_ZIP = os.path.join(ROOT, 'dist', 'full_archive', 'HAWRA_full_archive.zip')

class TestBundles(unittest.TestCase):
    def test_preprint_bundle_contents(self):
        self.assertTrue(os.path.exists(PREPRINT_ZIP))
        with zipfile.ZipFile(PREPRINT_ZIP, 'r') as z:
            names = set(z.namelist())
            required = {
                '00_docs/formalization/HAWRA-PQPE_formal_model.md',
                '05_data/results/hawra_plasmid.gb',
                '05_data/results/genetic_visualization/hawra_plasmid_map.png',
                '03_unified_simulator/results/spectral_sweep.csv',
                '01_genomics/raw_sequences/CDS/NOS.fasta',
            }
            for r in required:
                self.assertIn(r, names)
            for n in names:
                self.assertFalse(n.endswith('.pyc'))
                self.assertNotEqual(os.path.basename(n), '.DS_Store')

    def test_full_archive_bundle_contents(self):
        self.assertTrue(os.path.exists(FULL_ZIP))
        with zipfile.ZipFile(FULL_ZIP, 'r') as z:
            names = set(z.namelist())
            required = {
                'HAWRA/plasmid/hawra_plasmid.gb',
                'HAWRA/plasmid/hawra_plasmid_map.png',
                'HAWRA/protocol/SOP_ficus_regeneration.md',
            }
            for r in required:
                self.assertIn(r, names)
            for n in names:
                self.assertFalse(n.endswith('.pyc'))
                self.assertNotEqual(os.path.basename(n), '.DS_Store')

if __name__ == '__main__':
    unittest.main()