# Prevalent haemagglutinin 162–164 deletions in influenza B/Victoria virus influence fitness and virulence in vivo

## Overview
This repository contains the analysis code accompanying the study:

**"Prevalent haemagglutinin 162–164 deletions in influenza B/Victoria virus influence fitness and virulence in vivo"**

This study investigates how haemagglutinin (HA) deletions at residues 162–164 and associated amino acid substitutions influence the fitness, stability and virulence of influenza B/Victoria viruses. Influenza B viruses cause substantial respiratory disease and seasonal outbreaks. Following the COVID-19 pandemic, only the B/Victoria lineage persisted, yet the functional consequences of its emerging HA deletion variants remain unclear.

To address this gap, we combined phylodynamic analyses and reverse genetics to investigate the effects of HA deletions and substitutions across 19 B/Victoria recombinant viruses. Replication, receptor binding, HA stability and pathogenicity in mice were examined to understand the functional changes associated with the emergence and replacement of Clades V1A.1 and V1A.3.

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
