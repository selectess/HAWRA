    # HAWRA - Schéma de l'Architecture du projet HAWRA

Ce document présente le schéma de l'architecture HAWRA, en soulignant la hiérarchie de contrôle stricte où l'écosystème ARBOL est l'unique point d'entrée.
L'utilisateur interagit exclusivement via l'IDE ou le terminal ARBOL. Le code est compilé en un "Assembly Quantique", qui est ensuite interprété par le BioOS pour piloter le substrat biologique (PQPE).

```mermaid
graph TD
    subgraph HAWRA [HAWRA : Système de Calcul Bio-Intégré]

        subgraph ARBOL_Interface [ARBOL: Interface de Contrôle Unique]
            direction LR
            User([Utilisateur]) --> ARBOL_IDE(IDE / Terminal ARBOL)
            ARBOL_IDE -- écrit --> ArbolCode(Code ARBOL)
            ArbolCode --> Compiler(Bio-Compilateur)
            Compiler -- génère --> QuantumAssembly(Assembly Quantique)
        end

        subgraph BioOS_Layer [BioOS: Interface Bio-Computationnelle]
            QuantumAssembly -- est traduit par --> BioOS(Cœur du BioOS)
            BioOS -- Régule via signaux --> PQPE_Substrate
            Sensors -- Données brutes --> BioOS
        end

        subgraph PQPE_Substrate [PQPE: Substrat Biologique de Calcul]
            subgraph PQPE_Entity [Entité PQPE]
                Organism(Organisme Végétal Modifié)
                QuantumCircuits(Circuits Quantiques Biologiques)
                Metabolism(Métabolisme & Génétique)
            end
            Sensors(Capteurs: IR, Spectra...)
            Actuators(Actuateurs: LED, EM)
        end

        subgraph Data_Ecosystem [Écosystème de Données]
            Genomics(01_genomics: Données Génomiques)
            Results(Résultats d'expériences)
            Data(05_data: Stockage de données)
        end

    end

    %% Flux de Contrôle et de Données
    ARBOL_Interface -- 1. Transmet l'Assembly Quantique --> BioOS_Layer

    BioOS_Layer -- 2. Traduit en signaux --> Actuators
    Actuators -- 3. Stimule --> QuantumCircuits
    QuantumCircuits -- 4. État mesuré par --> Sensors
    Sensors -- 5. Remonte les données brutes --> BioOS_Layer

    BioOS_Layer -- 6. Traite et produit --> Results
    Results -- 7. Stocké dans --> Data
    Data -- 8. Visualisé dans --> ARBOL_Interface

    %% Liaisons de support
    Genomics -- Informe --> Metabolism
    Metabolism -- Support   e --> Organism
    Organism -- Héberge --> QuantumCircuits

    classDef software fill:#c9f,stroke:#333,stroke-width:2px;
    classDef hardware fill:#f9c,stroke:#333,stroke-width:2px;
    classDef data fill:#9cf,stroke:#333,stroke-width:2px;

    class BioOS_Layer, ARBOL_Interface software;
    class PQPE_Substrate, PQPE_Entity, Organism, QuantumCircuits, Metabolism, Sensors, Actuators hardware;
    class Data_Ecosystem, Genomics, Data, Results data;
```
