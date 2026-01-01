# HAWRA Project Structure

This document outlines the structure of the HAWRA project, providing an overview of the key directories and their contents.

## Root Directory

The root directory contains the following key files and directories:

- `README.md`: The main README file for the project.
- `STRUCTURE.md`: This file, which describes the project structure.
- `requirements.txt`: A list of the Python dependencies for the project.
- `scripts/`: A directory containing various scripts for analyzing the simulation results.
- `results/`: A directory containing the results of the simulations.

## Main Directories

The HAWRA project is organized into the following main directories:

- `00_docs/`: Contains documentation for the project, including papers, notes, and references.
- `01_genomics/`: Contains genomic data, including plasmid maps and raw sequences.
- `02_arbol_interface/`: Contains the source code for the Arbol IDE.
- `03_unified_simulator/`: Contains the source code for the main HAWRA simulator.
- `05_data/`: Contains various data files, including experimental logs and simulation results.
- `arbol/`: Contains the source code for the Arbol compiler.
- `bioos/`: Contains the source code for the BioOS.

## `03_unified_simulator` Directory

The `03_unified_simulator` directory is the heart of the HAWRA project. It contains the source code for the main simulator, as well as the following key files and directories:

- `main.py`: The main entry point for the simulator.
- `gene_control.bsim.json`: A sample BSIM script for running a gene control simulation.
- `src/`: The source code for the simulator, which is organized into the following subdirectories:
    - `hawra_simulator/`: The main simulator package.
        - `engines/`: The simulation engines for the biological system, the environment, and the quantum system.
        - `simulator.py`: The main simulator class.

## `scripts` Directory

The `scripts/` directory contains a number of scripts for analyzing the simulation results. Here are a few examples:

- `analyze_gene_control.py`: This script analyzes the results of the gene control simulation and generates a plot of the gene expression levels over time.
- `visualize_plasmid.py`: This script visualizes a plasmid map from a GenBank file.
- `visualize_results.py`: This script visualizes the results of a simulation.

## `results` Directory

The `results/` directory contains the results of the simulations. The results are typically saved in JSON format, and can be analyzed using the scripts in the `scripts/` directory.
