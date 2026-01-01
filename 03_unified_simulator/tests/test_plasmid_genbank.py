import os
import unittest
from Bio import SeqIO

class TestPlasmidGenBank(unittest.TestCase):
    def test_generated_genbank_integrity(self):
        out_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), '05_data', 'results', 'hawra_plasmid.gb')
        if not os.path.exists(out_path):
            import runpy
            runpy.run_path(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'scripts', 'create_genbank.py'), run_name='__main__')
        record = SeqIO.read(out_path, 'genbank')
        seq_str = str(record.seq)
        self.assertTrue('N' not in seq_str)
        expected_genes = ['psaA', 'CRY2', 'Luc', 'Lsi1', 'HSP70', 'PEPC']
        for g in expected_genes:
            has_prom = any(f.type == 'promoter' and f.qualifiers.get('label', [''])[0].endswith(g) for f in record.features)
            has_ins = any(f.type == 'misc_feature' and f.qualifiers.get('label', [''])[0].startswith('Insulator_') and f.qualifiers.get('label', [''])[0].endswith(g) for f in record.features)
            has_5utr = any(f.type == 'misc_feature' and f.qualifiers.get('label', [''])[0].startswith("5UTR_Kozak_") and f.qualifiers.get('label', [''])[0].endswith(g) for f in record.features)
            cds_list = [f for f in record.features if f.type == 'CDS' and f.qualifiers.get('gene', [''])[0] == g]
            has_term = any(f.type == 'terminator' and f.qualifiers.get('label', [''])[0].endswith(g) for f in record.features)
            self.assertTrue(has_prom)
            self.assertTrue(has_ins)
            self.assertTrue(has_5utr)
            self.assertTrue(len(cds_list) == 1)
            self.assertTrue(has_term)
            cds_seq = str(record.seq[int(cds_list[0].location.start):int(cds_list[0].location.end)])
            self.assertTrue(cds_seq.startswith('ATG'))

if __name__ == '__main__':
    unittest.main()