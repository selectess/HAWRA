# Architecture de BioOS et Interface Matérielle

Ce document formalise l'architecture du système BioOS, son interaction avec l'entité de traitement quantique phyto-synthétique (PQPE), et le contrat d'interface avec le matériel physique.

## 1. Vision d'Ensemble : Un Système Cyber-Physique

BioOS n'est pas un logiciel monolithique, mais un **système cyber-physique** composé de deux éléments principaux qui fonctionnent en tandem :

1.  **L'Orchestrateur BioOS (Partie Digitale)** : Un runtime logiciel, typiquement hébergé sur une plateforme embarquée comme un **NVIDIA Jetson**, qui sert de pont entre le monde digital et le monde biologique.
2.  **L'Exécuteur BioOS (Partie Biologique)** : La **PQPE** elle-même. C'est une entité biologique (plante) génétiquement modifiée pour héberger une logique computationnelle *in vivo* via un **plasmide computationnel**.

L'interaction entre ces deux composants forme une boucle de contrôle fermée : l'Orchestrateur envoie des stimuli physiques, et l'Exécuteur y répond en modifiant son état, qui est ensuite lu par l'Orchestrateur.

## 2. Le Flux d'Exécution : De l'Arbol à la Réponse Biologique

Le cycle de vie d'une instruction dans le système HAWRA est le suivant :

1.  **Programmation Haut Niveau** : Un scientifique ou un développeur écrit un programme en utilisant le langage **Arbol**, qui décrit une séquence d'opérations logiques.
2.  **Compilation** : Le compilateur Arbol transforme le script `.arbol` en un contrat bas niveau, le **BSIM** (`.bsim.json`). Ce fichier est un "bytecode" ou un "langage d'assemblage" pour le système HAWRA.
3.  **Orchestration (Traduction)** : L'**Orchestrateur BioOS** lit le contrat BSIM. Il ne l'exécute pas sur un CPU. À la place, il le **traduit** en une séquence de commandes pour les périphériques matériels auxquels il est connecté.
    *   **Exemple** : L'instruction BSIM `STIMULUS_APPLY 'RED_LIGHT' DURATION 2s` est traduite en un signal électrique qui active un driver de LED pour émettre une lumière rouge à 660nm pendant 2 secondes.
4.  **Stimulation Physique** : Les actionneurs (LEDs, lasers, bobines électromagnétiques) appliquent le stimulus physique à la PQPE. C'est le "billet d'entrée" de l'information dans le système biologique.
5.  **Exécution *In Vivo*** : La logique du **plasmide computationnel** au sein de la PQPE détecte le stimulus. Cela déclenche une cascade métabolique ou épigénétique prédéfinie, ce qui constitue l'"exécution" de l'instruction.
6.  **Réponse Biologique** : L'exécution *in vivo* produit un changement d'état mesurable dans la plante (ex: modification du potentiel électrique, émission de fluorescence, etc.).
7.  **Lecture (Sensing)** : Les capteurs connectés à l'Orchestrateur BioOS (caméras, électrodes, photomètres) détectent ce changement d'état.
8.  **Interprétation et Feedback** : L'Orchestrateur reçoit les données brutes des capteurs, les interprète, les enregistre, et peut potentiellement générer un nouveau rapport d'état (potentiellement sous forme de `.bsim` ou de logs) pour informer la suite du programme.

## 3. Le Contrat d'Interface Matérielle

L'Orchestrateur BioOS communique avec le matériel via une couche d'abstraction de drivers. Le contrat d'interface (API) pour ces drivers doit être clairement défini.

*   **Actionneurs (Sorties de l'Orchestrateur)** :
    *   `light.set(channel, intensity, wavelength)`
    *   `em_field.set(frequency, amplitude, waveform)`
    *   `fluidics.pump(volume, flow_rate)`
*   **Capteurs (Entrées vers l'Orchestrateur)** :
    *   `camera.read_frame()`
    *   `electrode.read_potential(channel)`
    *   `spectrometer.read_spectrum()`

## 4. Modes d'Exécution

BioOS supporte deux modes d'exécution :

*   **Mode Physique (Runtime)** : Le mode de fonctionnement normal, utilisant la boucle cyber-physique complète avec le Jetson, les drivers, et la PQPE.
*   **Mode Simulation (`fast_mode`)** : Un mode purement logiciel qui simule le comportement de la PQPE. Il consomme un contrat BSIM et produit des logs et des métriques simulés, sans aucune interaction matérielle. Ce mode est essentiel pour le développement, le débogage rapide et la validation de la logique des scripts Arbol.
