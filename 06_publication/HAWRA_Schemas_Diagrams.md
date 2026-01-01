# Schémas et Diagrammes du Projet HAWRA

Ces schémas sont conçus pour être intégrés dans le papier scientifique via Mermaid.js ou TikZ.

## 1. Architecture Globale ARBOL → Hardware
```mermaid
graph TD
    A[Script .arbol] -->|Compiler| B(Arbol Compiler Lark v0.3)
    B -->|Generate| C[BSIM JSON Contract]
    C -->|Execute| D{BioOS Kernel}
    D -->|Digital Twin| E[Unified Simulator]
    D -->|Physical Execution| F[Jetson Client API]
    F -->|Control| G[Hardware: LEDs/EM Field]
    G -->|Response| H[PQPE Biological Entity]
    H -->|Feedback| I[Sensors: Electrodes/LUC]
    I -->|Data| D
    E -->|Validation| D
```

## 2. Pipeline de Données Génomique
```mermaid
flowchart LR
    A[NCBI / Raw GB] -->|CDS Selection| B[Genetic Cassette Design]
    B -->|Annotation| C[HAWRA_FINAL_VALIDATED.gb]
    C -->|Simulation| D[Biological Engine]
    D -->|Result| E[P700 Coherence Curve]
```

## 3. Workflow d'Exécution BioOS
```mermaid
sequenceDiagram
    participant User
    participant BioOS
    participant Jetson
    participant Simulator
    User->>BioOS: Load .arbol
    BioOS->>BioOS: Compile to BSIM
    BioOS->>Jetson: Check Status (Port 5001)
    Jetson-->>BioOS: Ready
    BioOS->>Simulator: Run Parallel Digital Twin
    BioOS->>Jetson: Apply Stimulus (Light 660nm)
    Jetson->>BioOS: Return Sensor Data
    BioOS->>User: Display Result (WebApp)
```

## 4. Modèle de Couche d'Isolation (Silica Shield)
```mermaid
graph LR
    A[Thermal Noise] -->|Attenuation| B[Silica Deposition Lsi1]
    B -->|Protection| C[P700 Reaction Center]
    C -->|Output| D[Quantum State Maintenance]
```
