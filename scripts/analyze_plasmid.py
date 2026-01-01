from Bio import SeqIO

def analyze_genbank(file_path):
    """
    Analyzes a GenBank file and prints information about its features.
    """
    for record in SeqIO.parse(file_path, "genbank"):
        print(f"Description: {record.description}")
        print(f"Sequence Length: {len(record.seq)}")
        print("Features:")
        for feature in record.features:
            if feature.type != "source":
                print(f"  - Type: {feature.type}")
                print(f"    Location: {feature.location}")
                if "gene" in feature.qualifiers:
                    print(f"    Gene: {feature.qualifiers['gene'][0]}")
                if "product" in feature.qualifiers:
                    print(f"    Product: {feature.qualifiers['product'][0]}")

if __name__ == "__main__":
    genbank_file = "/Users/mehdiwhb/Desktop/HAWRA/01_genomics/plasmids/HAWRA_FINAL_VALIDATED.gb"
    analyze_genbank(genbank_file)
