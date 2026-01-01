import os
import unittest
from Bio import SeqIO

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
GB_PATH = os.path.join(ROOT, '05_data', 'results', 'hawra_plasmid.gb')
NOS_FASTA = os.path.join(ROOT, '01_genomics', 'raw_sequences', 'CDS', 'NOS.fasta')

class TestNOSExactness(unittest.TestCase):
    def test_nos_present_and_exact(self):
        self.assertTrue(os.path.exists(GB_PATH))
        self.assertTrue(os.path.exists(NOS_FASTA))
        nos_records = list(SeqIO.parse(NOS_FASTA, 'fasta'))
        self.assertTrue(len(nos_records) >= 1)
        nos_seq = str(nos_records[0].seq).upper()
        gb = SeqIO.read(GB_PATH, 'genbank')
        term_features = [f for f in gb.features if f.type == 'terminator']
        self.assertTrue(len(term_features) >= 1)
        for f in term_features:
            seq = str(f.extract(gb.seq)).upper()
            self.assertEqual(seq, nos_seq)

if __name__ == '__main__':
    unittest.main()