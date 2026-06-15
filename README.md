# Xenopus laevis gene expression explorer

Self-contained Shiny app bundle for browsing reanalysed Xenopus laevis
proteomics and RNA-seq data across developmental stages and adult tissues.

## What it shows

- **Mass spectrometry** (Van Itallie 2025 reanalysis): protein abundance
  across 11 developmental stages (~14,200 proteins, two TMT 11-plex
  replicates), plus keratin phosphosite trajectories from a Rep C
  re-search against the krt12.4.S-corrected database.
- **RNA-seq** (Session et al. 2016 reanalysis): TPM matrices for the
  developmental time course (11 stages: egg through st40) and adult
  tissue panel (13 tissues).
- **Combined view**: paired protein + RNA trajectories at the five
  developmental stages where both modalities sampled the same
  Nieuwkoop-Faber stage (Egg, St9, St12, St30, St40/41).

## Contents

```
xenopus_explorer_bundle/
├── app.R                                  # Shiny app (single file)
├── README.md                              # this file
├── abundance_gene_MD.tsv                  # FragPipe protein abundance
├── keratin_phospho_sites_perstage.tsv     # Rep C phosphosite summary
├── keratin_peptides_all_perstage.tsv      # Rep C peptide summary
├── session2016_dev_tpm.tsv                # RNA TPM, developmental
└── session2016_tissue_tpm.tsv             # RNA TPM, adult tissues
```

All data are pre-computed; the app does no real-time computation beyond
filtering, Z-score normalisation, and plotting.

## Run locally

```r
# From R, after cd-ing into this directory:
shiny::runApp(".")
# Or from anywhere:
shiny::runApp("path/to/xenopus_explorer_bundle")
```

Required R packages (all standard CRAN):

```r
install.packages(c("shiny", "dplyr", "tidyr", "ggplot2",
                   "plotly", "DT", "readr", "scales"))
```

## Resource requirements

- R 4.4 or newer.
- ~16 MB on disk (code + all data).
- ~500 MB - 1 GB RAM per active session.
- 1-2 CPU cores per session is sufficient.
- Light traffic (designed for a small lab + collaborators; typically
  1-3 concurrent users).

## Updating the data

Replace any of the five TSVs in this directory with newer versions; the
app picks them up on the next launch. Schemas to preserve:

- `abundance_gene_MD.tsv`: FragPipe `tmt-report/abundance_gene_MD.tsv`
  format (Index + metadata + 11 RepA stages + 11 RepB stages).
- `keratin_phospho_sites_perstage.tsv`,
  `keratin_peptides_all_perstage.tsv`: produced by
  `cluster_reanalysis/20_keratin_phospho_peptide_summary.py`.
- `session2016_dev_tpm.tsv`, `session2016_tissue_tpm.tsv`: produced by
  `session 2016/build_session2016_tpm.py`.
