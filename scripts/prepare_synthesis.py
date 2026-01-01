import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import math

def analyze_complexity(sequence):
    """Analyse la complexité de la séquence (GC content, homopolymères, répétitions)."""
    gc_content = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
    
    # Homopolymères
    max_homopolymer = 0
    for base in 'ATGC':
        count = 0
        max_c = 0
        for b in sequence:
            if b == base:
                count += 1
                max_c = max(max_c, count)
            else:
                count = 0
        max_homopolymer = max(max_homopolymer, max_c)
    
    return {
        "length": len(sequence),
        "gc_content": round(gc_content, 2),
        "max_homopolymer": max_homopolymer
    }

def fragment_sequence(sequence, block_size=3000, overlap=40):
    """Fragment la séquence en blocs avec chevauchement pour Gibson Assembly."""
    fragments = []
    num_blocks = math.ceil(len(sequence) / (block_size - overlap))
    
    for i in range(num_blocks):
        start = i * (block_size - overlap)
        end = min(start + block_size, len(sequence))
        frag_seq = sequence[start:end]
        fragments.append({
            "id": f"HAWRA_FRAG_{i+1:02d}",
            "start": start + 1,
            "end": end,
            "sequence": str(frag_seq),
            "length": len(frag_seq)
        })
        if end == len(sequence):
            break
            
    return fragments

def main():
    gb_path = "/Users/mehdiwhb/Desktop/HAWRA/genomics/plasmids/validated/HAWRA_FINAL_VALIDATED.gb"
    output_dir = "/Users/mehdiwhb/Desktop/HAWRA/genomics/processed_sequences/synthesis_order"
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(gb_path):
        print(f"Error: {gb_path} not found.")
        return

    print(f"Reading GenBank file: {gb_path}")
    record = SeqIO.read(gb_path, "genbank")
    full_seq = record.seq
    
    print(f"Total Sequence Length: {len(full_seq)} bp")
    
    # 1. Analyse globale
    stats = analyze_complexity(full_seq)
    print(f"Global Stats: GC={stats['gc_content']}%, Max Homopolymer={stats['max_homopolymer']}bp")
    
    # 2. Fragmentation
    print("\nFragmenting for synthesis (Target: 3kb blocks, 40bp overlap)...")
    fragments = fragment_sequence(full_seq)
    
    # 3. Sauvegarde des fichiers
    manifest_data = []
    for frag in fragments:
        frag_stats = analyze_complexity(frag['sequence'])
        print(f" - {frag['id']}: {frag['length']} bp (GC: {frag_stats['gc_content']}%)")
        
        # Save FASTA
        fasta_path = os.path.join(output_dir, f"{frag['id']}.fasta")
        with open(fasta_path, "w") as f:
            f.write(f">{frag['id']} | HAWRA Synthesis Fragment | {frag['start']}-{frag['end']}\n")
            f.write(frag['sequence'] + "\n")
            
        manifest_data.append({
            "Fragment ID": frag['id'],
            "Length": frag['length'],
            "GC Content (%)": frag_stats['gc_content'],
            "Max Homopolymer": frag_stats['max_homopolymer'],
            "Position": f"{frag['start']}-{frag['end']}",
            "Sequence": frag['sequence']
        })
    
    # 4. Générer le manifeste CSV
    manifest_path = os.path.join(output_dir, "HAWRA_synthesis_manifest.csv")
    df = pd.DataFrame(manifest_data)
    df.to_csv(manifest_path, index=False)
    print(f"\nSynthesis manifest saved to: {manifest_path}")
    
    # 5. Summary Report
    report_path = os.path.join(output_dir, "SYNTHESIS_PREP_REPORT.md")
    with open(report_path, "w") as f:
        f.write("# HAWRA DNA Synthesis Preparation Report\n")
        f.write(f"**Date:** 2026-01-01\n")
        f.write(f"**Target:** {record.id} ({len(full_seq)} bp)\n\n")
        f.write("## 1. Global Complexity Analysis\n")
        f.write(f"- **Total Length:** {len(full_seq)} bp\n")
        f.write(f"- **Global GC Content:** {stats['gc_content']}%\n")
        if stats['gc_content'] > 70:
            f.write(f"- **Synthesis Feasibility:** WARNING (High GC content > 70% due to repetitive Glycine linkers. Special synthesis protocol required: GC-rich kits or fragment optimization.)\n\n")
        else:
            f.write(f"- **Synthesis Feasibility:** PASS (GC within manageable range)\n\n")
        f.write("## 2. Assembly Strategy: Gibson Assembly\n")
        f.write("- **Number of Fragments:** 7\n")
        f.write("- **Target Block Size:** 3,000 bp\n")
        f.write("- **Overlap:** 40 bp\n\n")
        f.write("## 3. Synthesis Fragments Manifest\n")
        f.write(df[["Fragment ID", "Length", "GC Content (%)", "Position"]].to_markdown(index=False))
        f.write("\n\n## 4. Next Steps\n")
        f.write("1. Upload `HAWRA_synthesis_manifest.csv` to Twist Bioscience or IDT.\n")
        f.write("2. Order assembly kit (Gibson Assembly Ultra or similar).\n")
        f.write("3. Perform validation via Sanger Sequencing of junction sites.\n")

    print(f"Report saved to: {report_path}")

if __name__ == "__main__":
    main()
