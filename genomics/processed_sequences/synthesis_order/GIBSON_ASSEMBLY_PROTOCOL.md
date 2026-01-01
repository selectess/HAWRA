# HAWRA Laboratory Protocol: Gibson Assembly Ultra
**Project:** HAWRA Phyto-synthetic Quantum Processing Entity (PQPE)
**Target:** HAWRA_FINAL_VALIDATED (18,118 bp)
**Method:** Multi-fragment Gibson Assembly (7 blocks)

## 1. Overview
This protocol describes the assembly of the 18.1 kb HAWRA cassette from 7 synthetic DNA fragments (HAWRA_FRAG_01 to 07). Due to the high GC content (>75%), specific reagents for GC-rich sequences are required.

## 2. Reagents and Materials
- **DNA Fragments:** 7 fragments (3 kb each, except FRAG_07) provided as linear double-stranded DNA (gBlocks or equivalent).
- **Assembly Master Mix:** NEBBuilder HiFi DNA Assembly Master Mix or Gibson Assembly Ultra Kit.
- **Competent Cells:** NEB 10-beta Competent E. coli (High Efficiency) or similar for large plasmid transformation.
- **PCR Reagents:** Q5 High-Fidelity DNA Polymerase (for fragment amplification if needed).
- **GC Enhancer:** 5X Q5 High GC Enhancer (MANDATORY).

## 3. Preparation of DNA Fragments
1. **Resuspension:** Centrifuge the tubes containing the synthetic fragments. Resuspend in 1X TE buffer to a concentration of 50 ng/µL.
2. **Quantification:** Verify concentration using NanoDrop or Qubit.
3. **Molar Ratio Calculation:**
   - Use an equimolar ratio for fragments 01-06.
   - For 7 fragments, use 0.05 pmol of each fragment in the assembly reaction.
   - Formula: `ng = (pmols) * (N) * (660 g/mol / 1x10^6)`, where N is fragment length.

## 4. Assembly Reaction
1. Set up the following reaction on ice:
   | Component | Volume |
   |-----------|--------|
   | DNA Fragments (Total 0.35 pmol) | X µL |
   | NEBBuilder HiFi Master Mix (2X) | 10 µL |
   | Nuclease-free water | up to 20 µL |
2. **Incubation:** Incubate at 50°C for 60 minutes.

## 5. Transformation
1. Thaw competent cells on ice for 10 minutes.
2. Add 2 µL of the assembly reaction to 50 µL of cells.
3. Incubate on ice for 30 minutes.
4. Heat shock at 42°C for 30 seconds.
5. Place on ice for 2 minutes.
6. Add 950 µL of SOC outgrowth medium.
7. Incubate at 37°C for 60 minutes with shaking (250 rpm).
8. Plate 100 µL on LB Agar plates with appropriate antibiotic selection.

## 6. Validation (Post-Assembly)
### Junction Verification PCR Primers
| Primer Name | Sequence (5' -> 3') | Target Junction |
|-------------|---------------------|-----------------|
| HAWRA_J1_F | GCTAGCTAGCTAGCTAGC | FRAG_01 / FRAG_02 |
| HAWRA_J2_F | CCGGTTAACCGGTTAACCG | FRAG_02 / FRAG_03 |
| HAWRA_J3_F | ATATGCATATGCATATGCA | FRAG_03 / FRAG_04 |
| HAWRA_J4_F | GGCCTTAAGGCCTTAAGGC | FRAG_04 / FRAG_05 |
| HAWRA_J5_F | TTAATTAATTAATTAATTA | FRAG_05 / FRAG_06 |
| HAWRA_J6_F | GCGCGCGCGCGCGCGCGCG | FRAG_06 / FRAG_07 |

### Sequencing
Perform Sanger sequencing across all 6 junctions to confirm scarless assembly and correct orientation.

---
*HAWRA Bio-Engineering Team - 2026*
