# HAWRA - Hardware-Agnostic Wetware-Reliant Architecture

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Bio-Quantum](https://img.shields.io/badge/Bio--Quantum-Open%20Science-green)](docs/theory/HAWRA-PQPE_formal_model.md)

## 🌿 Overview
HAWRA is a pioneering open-source project dedicated to **Metabiotic Computing**: leveraging living biological substrates for ambient-temperature quantum processing. By integrating synthetic biology with quantum physics, HAWRA transforms the *Ficus elastica* into a Phyto-synthetic Quantum Processing Entity (PQPE).

## 🚀 Key Features
- **P700 Bio-Qubit**: Natural quantum coherence in Photosystem I centers.
- **Silica Shielding**: Genetic biomineralization for 300K decoherence protection.
- **Arbol DSL**: A high-level domain-specific language for bio-quantum programming.
- **Unified Simulator**: Multiphysics engine (Lindblad + Hill kinetics) for validation.

## 📁 Repository Structure
- `docs/`: Formal mathematical models, theory, and experimental protocols.
- `src/`: Core simulation engine and Arbol DSL interpreter.
- `genomics/`: Validated GenBank files for the PQPE plasmid.
- `examples/`: Reference simulation scenarios (.bsim) and test cases.
- `infrastructure/`: Docker and deployment configurations.

## 🛠 Quick Start
To reproduce the simulation results:

### Option A: Docker (Recommended)
1. **Build the environment**:
   ```bash
   docker build -t hawra-sim infrastructure/
   ```
2. **Run a simulation**:
   ```bash
   docker run hawra-sim python src/run.py --bsim examples/first_bloom.bsim.json --output results.json
   ```

### Option B: Native Setup
1. **Configure Environment**:
   ```bash
   conda activate hawra_theory
   export PYTHONPATH=$PYTHONPATH:$(pwd)/src
   ```
2. **Run a simulation**:
   ```bash
   python3 src/run.py --bsim examples/first_bloom.bsim.json --output results.json
   ```

## 📖 Documentation
- [Formal Model](docs/theory/HAWRA-PQPE_formal_model.md)
- [Genetic Engineering Protocol](docs/protocols/SOP_ficus_regeneration.md)
- [Scientific Foundations](docs/theory/PQPE_Math_Bio_Engineering.md)

## 🤝 Contributing
Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## 📜 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🏷 Citation
If you use this work in your research, please cite it using:
```bibtex
@software{hawra_2026,
  author = {HAWRA Team},
  title = {HAWRA: Metabiotic Quantum Computing Architecture},
  year = {2026},
  url = {https://github.com/your-repo/hawra}
}
```
