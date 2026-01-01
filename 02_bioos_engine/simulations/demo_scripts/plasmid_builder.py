from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Charger les séquences depuis les fichiers locaux
seqs = {}
local_files = {
    "CRY2": "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/raw_sequences/CRY2_NM_100320.4.fasta",
    "HSP70": "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/raw_sequences/HSP70_NM_112093.3.fasta",
    "LUC": "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/raw_sequences/LUC_M15077.1.fasta",
    "Lsi1": "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/raw_sequences/Lsi1_SIT1_AB222272.1.fasta",
    "PEPC": "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/raw_sequences/PEPC1_X13660.1.fasta",
    "psaA": "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/raw_sequences/psaA_NC_000932.1.fasta"
}

for name, path in local_files.items():
    seqs[name] = SeqIO.read(path, "fasta")

# Mutation S284T pour LUC_red
luc_red_seq = seqs["LUC"].seq
luc_red_seq = luc_red_seq[:850] + Seq("ACT") + luc_red_seq[853:]

# Assemblage simplifié
cassette = (
    seqs["psaA"].seq[:100] +  # Promoteur minimal
    seqs["CRY2"].seq + 
    luc_red_seq + 
    seqs["Lsi1"].seq
)

final = SeqRecord(cassette, id="HAWRA_LOCAL_BUILD", description="Assemblage avec séquences locales - Test")
final.annotations["molecule_type"] = "DNA"
SeqIO.write(final, "HAWRA_PLASMID_LOCAL_TEST.gb", "genbank")

print(f"Plasmide assemblé localement : {len(cassette)} bp")
