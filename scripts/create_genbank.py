from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

# 1. Définir la séquence de base (pGreenII0229 + cassette)
backbone_seq = ("ATGC" * 2500)

# Séquences fictives pour les gènes de la cassette HAWRA
def read_fasta_sequence(path):
    try:
        records = list(SeqIO.parse(path, "fasta"))
        if records:
            seq = str(records[0].seq)
            return seq if len(seq) > 0 else None
    except Exception:
        pass
    return None

def prefer_first_valid(paths):
    for p in paths:
        s = read_fasta_sequence(p)
        if s:
            return s
    return None

def read_optional_3utr(gene_name):
    candidates = [
        f"01_genomics/raw_sequences/CDS/3UTR/{gene_name}.fasta",
        f"01_genomics/raw_sequences/CDS/{gene_name}_3UTR.fasta",
    ]
    return prefer_first_valid(candidates)

psaa_seq = prefer_first_valid([
    "01_genomics/raw_sequences/CDS/psaA_CDS_CORRECTED.fasta",
    "01_genomics/raw_sequences/CDS/psaA_CDS.fasta",
    "01_genomics/raw_sequences/CDS/psaA.fasta",
]) or ("ATG" + "N" * 1497)
cry2_seq = read_fasta_sequence("01_genomics/raw_sequences/CDS/CRY2_CDS.fasta") or ("ATG" + "N" * 1497)
luc_seq = read_fasta_sequence("01_genomics/raw_sequences/CDS/LUC_CDS.fasta") or ("ATG" + "N" * 1497)
lsi1_seq = read_fasta_sequence("01_genomics/raw_sequences/CDS/Lsi1_SIT1_CDS.fasta") or ("ATG" + "N" * 1497)
hsp70_seq = read_fasta_sequence("01_genomics/raw_sequences/CDS/HSP70_CDS.fasta") or ("ATG" + "N" * 1497)
pepc_seq = read_fasta_sequence("01_genomics/raw_sequences/CDS/PEPC1_CDS.fasta") or ("ATG" + "N" * 1497)
promoter_seq = prefer_first_valid([
    "01_genomics/raw_sequences/CDS/35S.fasta",
    "01_genomics/raw_sequences/35S.fasta"
]) or "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC"
insulator_seq = "AT" * 38
kozak_5utr_seq = "ACC"
terminator_seq = prefer_first_valid([
    "01_genomics/raw_sequences/CDS/NOS.fasta",
    "01_genomics/raw_sequences/NOS.fasta"
]) or ("AT" * 60)

# Assembler la cassette
cassette_parts = []
for cds in [
    ("psaA", psaa_seq),
    ("CRY2", cry2_seq),
    ("Luc", luc_seq),
    ("Lsi1", lsi1_seq),
    ("HSP70", hsp70_seq),
    ("PEPC", pepc_seq),
]:
    cassette_parts.append(promoter_seq)
    cassette_parts.append(insulator_seq)
    cassette_parts.append(kozak_5utr_seq)
    cassette_parts.append(cds[1])
    cassette_parts.append(terminator_seq)
cassette_seq = "".join(cassette_parts)

# Insérer la cassette dans le backbone (à une position arbitraire)
insertion_point = 1000
full_seq = backbone_seq[:insertion_point] + cassette_seq + backbone_seq[insertion_point:]
plasmid_seq = Seq(full_seq)

# 2. Créer l'objet SeqRecord
record = SeqRecord(
    plasmid_seq,
    id="pHAWRA01",
    name="pHAWRA",
    description="Plasmide computationnel pGreenII0229 optimisé pour exécution HAWRA",
    annotations={"molecule_type": "DNA", "topology": "circular"},
)

# 3. Définir les features (gènes de la cassette)
features = []
cursor = insertion_point
for gene_name, gene_seq in [
    ("psaA", psaa_seq),
    ("CRY2", cry2_seq),
    ("Luc", luc_seq),
    ("Lsi1", lsi1_seq),
    ("HSP70", hsp70_seq),
    ("PEPC", pepc_seq),
]:
    features.append(SeqFeature(FeatureLocation(cursor, cursor + len(promoter_seq), strand=1), type="promoter", qualifiers={"label": f"P35S_{gene_name}", "note": "CaMV 35S promoter"}))
    cursor += len(promoter_seq)
    features.append(SeqFeature(FeatureLocation(cursor, cursor + len(insulator_seq), strand=1), type="misc_feature", qualifiers={"label": f"Insulator_{gene_name}", "note": "RiboJ/insulator"}))
    cursor += len(insulator_seq)
    features.append(SeqFeature(FeatureLocation(cursor, cursor + len(kozak_5utr_seq), strand=1), type="misc_feature", qualifiers={"label": f"5UTR_Kozak_{gene_name}", "note": "Plant Kozak 5' UTR"}))
    cursor += len(kozak_5utr_seq)
    features.append(SeqFeature(FeatureLocation(cursor, cursor + len(gene_seq), strand=1), type="CDS", qualifiers={"gene": gene_name, "product": f"{gene_name}_protein"}))
    cursor += len(gene_seq)
    opt_3utr = read_optional_3utr(gene_name)
    if opt_3utr:
        features.append(SeqFeature(FeatureLocation(cursor, cursor + len(opt_3utr), strand=1), type="misc_feature", qualifiers={"label": f"3UTR_{gene_name}", "note": "Plant 3' UTR (optional)"}))
        cursor += len(opt_3utr)
    features.append(SeqFeature(FeatureLocation(cursor, cursor + len(terminator_seq), strand=1), type="terminator", qualifiers={"label": f"NOS_{gene_name}", "note": "NOS terminator"}))
    cursor += len(terminator_seq)
record.features = features

# 4. Écrire le fichier GenBank
output_file = "05_data/results/hawra_plasmid.gb"
with open(output_file, "w") as f:
    SeqIO.write(record, f, "genbank")

print(f"Fichier GenBank du plasmide HAWRA sauvegardé dans {output_file}")
