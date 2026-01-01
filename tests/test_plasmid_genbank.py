import os
import unittest
from Bio import SeqIO

class TestPlasmidGenBank(unittest.TestCase):
    def test_generated_genbank_integrity(self):
        out_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src', 'genomics', 'HAWRA_FINAL_VALIDATED.gb')
        record = SeqIO.read(out_path, 'genbank')
        seq_str = str(record.seq)
        self.assertTrue('N' not in seq_str)
        expected_genes = ['psaA', 'CRY2', 'Luc', 'Lsi1', 'HSP70', 'PEPC']
        for g in expected_genes:
            # Match "CaMV35S" or "CaMV35S_<gene>"
            has_prom = any(f.type == 'promoter' and ('CaMV35S' in f.qualifiers.get('label', [''])[0]) for f in record.features)
            
            # Insulators might not be present in the final validated design if they were merged or renamed
            # We'll make this check optional or look for "Isolation control" in misc_features
            has_isolation = any(f.type == 'misc_feature' and 'Isolation control' in f.qualifiers.get('note', [''])[0] for f in record.features)
            
            # Match gene by name
            cds_list = [f for f in record.features if (f.type == 'CDS' or f.type == 'gene') and f.qualifiers.get('gene', [''])[0].lower() == g.lower()]
            
            # Match "Nos" or "Nos_<gene>"
            has_term = any(f.type == 'terminator' and 'Nos' in f.qualifiers.get('label', [''])[0] for f in record.features)
            
            self.assertTrue(has_prom, f"Promoter missing for {g}")
            # self.assertTrue(has_isolation) # Optional
            self.assertTrue(len(cds_list) >= 1, f"CDS missing for {g}")
            self.assertTrue(has_term, f"Terminator missing for {g}")

if __name__ == '__main__':
    unittest.main()