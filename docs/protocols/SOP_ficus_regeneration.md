# Standard Operating Procedure (SOP): HAWRA-Ficus-G0 Generation
# Version 1.0

## 1. Objectif
Ce document fournit des instructions détaillées pour la génération d'une plante *Ficus elastica* transgénique (HAWRA-Ficus-G0) en utilisant la transformation médiée par *Agrobacterium tumefaciens*, sur la base du protocole théorique `ficus_regeneration_protocol.md`.

## 2. Sécurité
- Manipuler tous les organismes génétiquement modifiés (OGM) conformément aux réglementations locales et aux bonnes pratiques de laboratoire.
- Utiliser un équipement de protection individuelle (EPI) approprié : blouse, gants, lunettes de sécurité.
- Toutes les manipulations d'*Agrobacterium* et de tissus végétaux doivent être effectuées dans une hotte à flux laminaire stérile.

---

## 3. Matériels et Réactifs

### 3.1. Souches Bactériennes
- *Agrobacterium tumefaciens*, souche EHA105 (ou LBA4404).
- *E. coli*, souche DH5α (pour le clonage plasmidique).

### 3.2. Vecteurs
- Vecteur binaire pour Agrobacterium, série pCAMBIA (ex: pCAMBIA1301).
- Plasmide `HAWRA_FINAL_VALIDATED.gb` (séquence à synthétiser).

### 3.3. Milieux de Culture
- Milieu LB (Luria-Bertani) pour *E. coli* et *Agrobacterium*.
- Milieu MS (Murashige & Skoog), sels et vitamines.
- Milieu de co-culture (MS + 20 g/L sucrose, 100 µM acetosyringone).
- Milieu de sélection (MS + antibiotiques).
- Milieu de pousse (MS + BAP, Kinetin).
- Milieu d'enracinement (MS + IBA).

### 3.4. Antibiotiques et Produits Chimiques
- Kanamycine (pour sélection *E. coli* et plante).
- Rifampicine (pour sélection *Agrobacterium*).
- Céfotaxime (pour éliminer *Agrobacterium* post-infection).
- Acetosyringone.
- Hormones végétales : BAP (6-Benzylaminopurine), Kinetin, IBA (Indole-3-butyric acid).
- Agar, Sucrose, MES hydrate.

---

## 4. Procédure Détaillée

### Phase I : Préparation du Vecteur (Durée estimée : 2 semaines)

**Étape 1.1 : Synthèse et Clonage du Circuit HAWRA**
1.  **Commander la synthèse** de la séquence d'ADN correspondant à `HAWRA_FINAL_VALIDATED.gb`, avec codons optimisés pour *Ficus elastica*. Inclure des sites de restriction compatibles avec le vecteur pCAMBIA.
2.  **Digérer** le vecteur pCAMBIA et l'insert HAWRA avec les enzymes de restriction appropriées.
3.  **Ligaturer** l'insert dans le vecteur pCAMBIA.
4.  **Transformer** des *E. coli* DH5α compétents avec le produit de ligature.
5.  **Plaquer** sur un milieu LB + Kanamycine. Incuber à 37°C pendant une nuit.
6.  **Sélectionner** plusieurs colonies, les cultiver en milieu liquide LB + Kanamycine.
7.  **Extraire** l'ADN plasmidique (miniprep).
8.  **Vérifier** l'insertion correcte par digestion de restriction et séquençage Sanger.

**Étape 1.2 : Transformation d'Agrobacterium**
1.  Préparer des cellules d'*Agrobacterium* EHA105 électrocompétentes.
2.  **Transformer** les cellules avec le plasmide pCAMBIA-HAWRA vérifié.
3.  **Plaquer** sur un milieu LB + Kanamycine + Rifampicine. Incuber à 28°C pendant 2-3 jours.
4.  **Confirmer** la présence du plasmide dans les colonies d'*Agrobacterium* par PCR.
5.  Préparer des stocks glycérols de la souche validée et les conserver à -80°C.

---
*(Les phases suivantes seront détaillées ultérieurement)*

---

## 5. Variante : Protocole d'Optimisation

En cas de faible efficacité de transformation (échecs répétés à l'étape 2), une variante optimisée du protocole peut être utilisée. Cette variante vise à augmenter l'efficacité du transfert d'ADN par *Agrobacterium*.

- **Modification Clé :** Augmenter la concentration d'**acetosyringone** dans le milieu de co-culture de 100 µM à **200 µM**.
- **Justification :** Une concentration plus élevée d'acetosyringone peut induire plus fortement les gènes de virulence (vir) d'*Agrobacterium*, améliorant ainsi le transfert de l'ADN-T vers la cellule végétale.

Cette modification n'affecte que le milieu de co-culture et les étapes immédiatement associées.

---

## 6. Simulations Procédurales (Visualisation)

Pour aider à visualiser les défis et les points de défaillance potentiels de ce protocole, deux simulations procédurales sont disponibles. Elles suivent le parcours d'un **seul explant** à travers chaque étape critique.

### 6.1. Simulation du Protocole Standard

Pour lancer la simulation du protocole standard (100 µM acetosyringone), exécutez la commande suivante :

```bash
python3 03_bioos/simulations/sop_procedural_simulation.py
```

### 6.2. Simulation du Protocole Optimisé

Pour visualiser l'impact de l'optimisation, lancez la simulation correspondante (200 µM acetosyringone) :

```bash
python3 03_bioos/simulations/sop_procedural_simulation_optimized.py
```


