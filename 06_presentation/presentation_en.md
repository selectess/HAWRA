# HAWRA - Final Presentation of the HAWRA Project

## 1. Introduction

The HAWRA project aims to develop a biological computing platform capable of executing complex algorithms by leveraging the principles of synthetic biology and quantum mechanics. This document presents the model architecture, simulation results, and the orchestration capabilities developed.

**HAWRA** is an integrated system consisting of three pillars:
1.  **PQPE (Phyto-synthetic Quantum Processing Entity)**: The biological hardware (wetware).
2.  **BioOS**: The operating system that orchestrates it.
3.  **ARBOL System**: The programming ecosystem for BioOS.

## 2. Model Architecture

The multiphysics simulator is at the heart of the HAWRA project. It integrates three specialized engines to model the system's complex behavior:

*   **Environmental Engine**: Simulates external conditions, particularly light pulses that serve as stimuli for the system.
*   **Biological Engine**: Models the cell's genetic response, specifically the P700 protein concentration in response to light.
*   **Quantum Engine**: Describes the quantum state of P700, acting as a qubit, and its decoherence, which results in the emission of luciferase signals.

## 3. Simulation Results

### 3.1. Response to a Single Light Pulse

A initial simulation was performed to validate the system's response to a series of light pulses. The graph below shows the dynamics of P700 concentration and luciferase signals in response to light stimuli.

![Pulse Response](/static/images/simulation.png)

### 3.2. Parameter Sweep

To explore the model's robustness and sensitivity, a parameter sweep was conducted. Simulations were repeated for different values of P700 protein synthesis and degradation rates. The results show how these parameters influence the system dynamics.

#### Low Degradation, Low Synthesis (`deg_0.05_syn_0.2`)

![Low Degradation, Low Synthesis](/static/images/deg_0.05_syn_0.2.png)

In this regime, P700 concentration increases slowly and reaches a moderate plateau. The system is stable but not very reactive.

#### Medium Degradation, High Synthesis (`deg_0.1_syn_1.0`)

![Medium Degradation, High Synthesis](/static/images/deg_0.1_syn_1.0.png)

With high synthesis, P700 concentration increases rapidly and reaches high levels. The system is very reactive but may be subject to saturation.

#### High Degradation, Medium Synthesis (`deg_0.2_syn_0.5`)

![High Degradation, Medium Synthesis](/static/images/deg_0.2_syn_0.5.png)

High degradation leads to a rapid decrease in P700 concentration after each pulse. The system is highly dynamic and quickly resets its state, which can be desirable for fast computations.

## 4. ARBOL Orchestration

To automate and manage complex simulation campaigns, orchestration scripts were developed using the ARBOL framework. These scripts allow for:

*   **Defining and executing single experiments** with specific configurations (`run_pulse_experiment.py`).
*   **Automating parameter sweeps** to explore the solution space (`run_parameter_sweep.py`).

This approach facilitates reproducible research and large-scale analysis of model behavior.

## 5. Conclusion and Future Outlook

The HAWRA project has developed a functional multiphysics simulator and robust orchestration tools. Simulation results validate the model architecture and demonstrate its ability to reproduce complex biological and quantum behaviors.

Next steps could include:

*   Integrating more detailed models for biological and quantum processes.
*   Developing machine learning algorithms to analyze simulation results.
*   Experimental validation of the model in the laboratory.
