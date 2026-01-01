import os
from Bio import Entrez, SeqIO

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT = os.path.join(ROOT, '01_genomics', 'raw_sequences', 'CDS', 'NOS.fasta')

def fetch_pBI121_nos(accession='AF485783'):
    Entrez.email = 'example@example.com'
    handle = Entrez.efetch(db='nuccore', id=accession, rettype='gb', retmode='text')
    record = SeqIO.read(handle, 'genbank')
    handle.close()
    # heuristique: prendre la feature 'regulatory' avec /regulatory_class="terminator" et label/note mentionnant 'NOS'
    for f in record.features:
        if f.type == 'regulatory':
            rc = f.qualifiers.get('regulatory_class', [''])[0].lower()
            lbl = f.qualifiers.get('label', [''])[0].lower()
            note = f.qualifiers.get('note', [''])[0].lower()
            if rc == 'terminator' and ('nos' in lbl or 'nos' in note):
                seq = str(f.extract(record.seq))
                return seq
    # fallback: chercher terminator annoté 'NOS terminator' générique
    for f in record.features:
        if f.type == 'terminator':
            lbl = f.qualifiers.get('label', [''])[0].lower()
            note = f.qualifiers.get('note', [''])[0].lower()
            if 'nos' in lbl or 'nos' in note:
                seq = str(f.extract(record.seq))
                return seq
    return None

def main():
    seq = fetch_pBI121_nos()
    if not seq:
        raise SystemExit('NOS terminator introuvable depuis NCBI AF485783')
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as f:
        f.write('>NOS_terminator_pBI121\n')
        # wrap fasta 80 cols
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')
    print('NOS.fasta sauvegardé:', OUT)

if __name__ == '__main__':
    main()