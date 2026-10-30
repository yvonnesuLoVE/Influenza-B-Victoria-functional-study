# Prevalent haemagglutinin 162–164 deletions in influenza B/Victoria virus influence fitness and virulence in vivo

## Overview
This repository contains the analysis code accompanying the study:

**"Prevalent haemagglutinin 162–164 deletions in influenza B/Victoria virus influence fitness and virulence in vivo"**

This study investigates how haemagglutinin (HA) deletions at residues 162–164 and associated amino acid substitutions influence the fitness, stability and virulence of influenza B/Victoria viruses. Influenza B viruses cause substantial respiratory disease and seasonal outbreaks. Following the COVID-19 pandemic, only the B/Victoria lineage persisted, yet the functional consequences of its emerging HA deletion variants remain unclear.

To address this gap, we combined phylodynamic analyses and reverse genetics to investigate the effects of HA deletions and substitutions across 19 B/Victoria recombinant viruses. Replication, receptor binding, HA stability and pathogenicity in mice were examined to understand the functional changes associated with the emergence and replacement of Clades V1A.1 and V1A.3.


## Key Findings

- Clade V1A.1 viruses carrying a two-amino acid deletion and the substitutions D129G, I180V and R498K exhibited enhanced replication and virulence compared with the ancestral B/Brisbane/60/2008 virus lacking deletions.
- V1A.1 mutants carrying the dual combinations 129D+180I, 2aa reinsertion+129D, or 2aa reinsertion+498R caused greater mortality in mice than the corresponding single mutants.
- The prevailing Clade V1A.3, carrying a three-amino acid deletion and the K136E substitution, displaced V1A.1 and showed increased pathogenicity in mice.
- Reinsertion of deleted residues reduced overall HA stability, whereas combinations of two- or three-amino acid deletions with specific HA substitutions restored or enhanced stability.
- Across 19 B/Victoria recombinant viruses, greater thermal HA stability was strongly associated with higher lung viral loads, weight loss and mortality, identifying thermal stability as a key correlate of pathogenicity.

---

## Data Sources

- Influenza B/Victoria sequence data used for phylodynamic analyses of HA deletion variants and the emergence and replacement of viral clades.
- Experimental measurements from 19 B/Victoria recombinant viruses generated using reverse genetics, including viruses with HA deletions, residue reinsertions and associated substitutions.
- Viral replication, receptor-binding, pH stability and thermal stability measurements used to characterize the functional effects of HA variation.
- Mouse infection measurements, including lung viral loads, weight loss and mortality, used to assess pathogenicity and its association with thermal HA stability.

The experimental workbooks, sequencing reads and reference sequences are not bundled with the code package. Dataset accession details and availability information should accompany the study release.

---

## Repository Contents

- `FluB_public_code/scripts/assay_analysis.R`: combined analysis script for pH stability, thermal stability and receptor binding. It includes data preparation, summary statistics, heatmaps, linear mixed models, Tukey-adjusted pairwise comparisons and receptor-binding plots.
- `FluB_public_code/scripts/sequencing.sh`: staged sequencing workflow covering read trimming, reference indexing, mapping, consensus generation, quality control, optional taxonomic classification and consensus FASTA header formatting.
- `FluB_public_code/scripts/plot_tree.R`: optional utility for visualizing supplied phylogenetic trees with tip metadata. It does not perform phylogenetic inference or phylodynamic analyses.
- `FluB_public_code/config/`: configurable input and output locations, analysis settings and sequencing tool settings.
- `FluB_public_code/examples/`: example sample identifiers, reference identifiers, virus metadata and heatmap panel definitions.
- `FluB_public_code/README.md`: software requirements, input formats and instructions for running the scripts.
- `FluB_public_code/METHODS_AND_CHANGES.md` and `FluB_public_code/VALIDATION.md`: analysis processing choices, script consolidation details and validation scope.

---

## Implications

The findings demonstrate that combinations of HA deletions and substitutions can confer fitness advantages to emerging influenza B/Victoria variants, contributing to global selective sweeps, evolutionary success and the lineage's long-term persistence. The strong association between thermal HA stability and pathogenicity highlights the importance of understanding how HA variation shapes viral behaviour, with implications for vaccine efficacy and influenza control.

## Citation
If you use this repository in your research, please cite:

**[Your Full Citation Here]**
