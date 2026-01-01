# HAWRA - Hardware-Agnostic Wetware-Reliant Architecture

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Bio-Quantum](https://img.shields.io/badge/Bio--Quantum-Open%20Science-green)](00_docs/scientific/HAWRA_Scientific_Report_Full.md)
[![Preprint](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17908061-blue)](https://doi.org/10.5281/zenodo.17908061)
[![Validation](https://img.shields.io/badge/Status-VALIDATED%20%E2%9C%85-brightgreen)](04_validation_benchmarks/CERTIFICATION_REPORT.md)
[![GitHub Pages](https://img.shields.io/badge/Live-Documentation-orange)](https://selectess.github.io/HAWRA/)

## 🌿 Overview
HAWRA is the first pioneering open-source project dedicated to **Metabiotic Computing**: leveraging living biological substrates for ambient-temperature quantum processing. By integrating synthetic biology with quantum physics, HAWRA transforms the *Ficus elastica* into a **Phyto-synthetic Quantum Processing Entity (PQPE)**.

**"On ne rêve plus, on compile."**

---

## 🏆 Current Status: LAB READY 🧪
The HAWRA architecture has successfully passed the end-to-end numerical validation for advanced quantum benchmarks and is now ready for biological synthesis.
- **Quantum Algorithms Validated**: 
  - **First Bloom**: Initial quantum state preparation.
  - **Grover's Search**: Amplitude amplification in biological substrates.
  - **Deutsch-Jozsa**: Determination of biological function parity.
- **Synthesis Ready**: Full fragmentation (7 blocks) and [Gibson Assembly Protocol](genomics/processed_sequences/synthesis_order/GIBSON_ASSEMBLY_PROTOCOL.md) finalized for the 18.1 kb cassette.
- **Full Report**: [View the Latest Certification Report](04_validation_benchmarks/CERTIFICATION_REPORT.md)

------

## 🌐 Project Presentation & Live Demo
Explore our professional web interface for technical documentation and interactive simulation:
👉 **[https://selectess.github.io/HAWRA/](https://selectess.github.io/HAWRA/)**

---

## 📁 Repository Structure (Academic Standard)
The repository is organized following the 00-04 academic convention for clarity and reproducibility:

- **[00_docs/](00_docs/)**: Full scientific reports, theoretical models (PQPE), and Roadmap 2026.
  - 📄 [Scientific Report (PDF/MD)](00_docs/scientific/HAWRA_Scientific_Report_Full.md)
- **[01_genomics/](01_genomics/)**: Validated GenBank files and CRISPR-Cas9 protocols.
  - 🧬 [pHAWRA_v5.8_Validated.gb](01_genomics/plasmids/validated/HAWRA_FINAL_VALIDATED.gb)
- **[02_bioos_engine/](02_bioos_engine/)**: Core Biological Operating System, Arbol DSL compiler, and multiphysics engine.
  - ⚙️ [Arbol Compiler](02_bioos_engine/arbol_compiler/)
- **[03_quantum_simulation/](03_quantum_simulation/)**: Lindblad decoherence models and Bloch sphere dynamics simulations.
  - ⚛️ [Lindblad Solver](03_quantum_simulation/lindblad_master_equation.py)
- **[04_validation_benchmarks/](04_validation_benchmarks/)**: "First Bloom" simulation scenarios and benchmark results.
  - 📊 [Benchmark Scenario](04_validation_benchmarks/first_bloom.arbol)

---

## 🔬 Scientific Validation
The project architecture has been numerically validated using a **Digital Twin** approach, coupling quantum master equations with biological metabolic flux analysis.

### Theoretical Pillars:
1. **P700 Coherence**: T2 prolongation via biomineralized Silica cage.
2. **Lindblad-Hill Coupling**: Unified mathematical framework for bio-quantum state evolution.
3. **BSIM Architecture**: Biological Instruction Set for molecular orchestration.

---

## 🚀 Quick Start (Reproduction)

### Option A: Docker (Recommended)
1. **Build the environment**:
   ```bash
   docker build -t hawra-sim infrastructure/
   ```
2. **Run the "First Bloom" simulation**:
   ```bash
   docker run hawra-sim python 02_bioos_engine/simulations/multiphysics_simulator/validate_simulation.py
   ```

### Option B: Native Setup (Conda)
1. **Configure Environment**:
   ```bash
   conda activate hawra_theory
   export PYTHONPATH=$PYTHONPATH:$(pwd)/02_bioos_engine
   ```
2. **Run a simulation**:
   ```bash
   python3 02_bioos_engine/simulations/multiphysics_simulator/validate_simulation.py
   ```

---

## 📖 Key Technical Concepts
- **P700 Bio-Qubit**: Native quantum register in Photosystem I reaction centers.
- **Silica Shielding**: Genetic engineering of the Lsi1 transporter for T2 coherence prolongation at 300K.
- **Arbol DSL**: Domain-specific language for programming biological quantum logic.
- **Lindblad + Hill Coupling**: Unified modeling of quantum decoherence and enzymatic reaction kinetics.

---

## 🤝 Contributing & Open Science
This project follows **Open Science** principles. We invite researchers from quantum physics, synthetic biology, and computer science to contribute.
Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct.

---

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🏷 Citation
If you use HAWRA in your research, please cite:
```bibtex
@preprint{wahbi_hawra_2025,
  author = {Wahbi, Mehdi},
  title = {HAWRA: First Plant-Based Quantum OS with Native Machine Learning},
  year = {2025},
  doi = {10.5281/zenodo.17908061},
  url = {https://github.com/selectess/HAWRA}
}
```
