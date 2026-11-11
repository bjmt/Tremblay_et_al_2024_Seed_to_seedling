# `figures/`

The scripts in this folder generate the plots shown in the main and supplementary figures. They are run after the scripts in `analyses/`, from the same working directory, and each script writes its plots to a sub-folder of the same name (e.g. `Figure4/`). The inputs of each script are listed at the top. Most are outputs of the `preprocess/` and `analyses/` scripts; the remaining ones are the RT-qPCR, ChIP-qPCR and luciferase measurements in `data/qPCR/`, `data/ChIP-qPCR/` and `data/luciferase/`, `Functions.R` (from `analyses/`), `genes.bed` in `data/annotations/`, the CPC2 results in `data/RNAseq-transcripts/` (obtained with the [CPC2 web server](https://cpc2.gao-lab.org) from the FASTA files written by `Figure4.R`), and the external files listed in the main README.

The genome browser tracks were made using the `trackvizir` package included in this folder, which can be installed with `R CMD INSTALL trackvizir`.

| Script | Panels |
|---|---|
| `Figure1.R` | 1b, 1d, S1a-c, S1f-r, S11, S12 |
| `Figure2.R` | 2a-e, S2a-f |
| `Figure3.R` | 3a-d, S3a, S3b, S3d-g |
| `Figure4.R` | 4a-e, S4a-j, S13a-h |
| `Figure5.R` | 5a, 5b, 5d-g, S5a, S5c-j, S6a-k |
| `Figure6.R` | 6a-d, 6g, 6h, S7a-i, S8b, S8e, S8f |
| `Figure7.R` | 7a, 7f, S9a-i, S9k, S9l, S10b, S10d, S10f, S10g, S10j, S10k |
| `SupplementaryFigure10hi.R`* | S10h, S10i (table contents) |
| `SupplementaryFigure14.R`* | S14a-d |
| `SupplementaryFigure15.R`* | S15a, S15b |
| `SupplementaryFigure16.R`* | S16 |

\* The original code for these panels was lost; these scripts recreate them to reproduce the published panels near exactly.
