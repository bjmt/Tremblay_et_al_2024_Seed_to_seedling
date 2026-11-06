# `analyses/`

The scripts in this folder reproduce the main analyses from the manuscript, and generate the processed data used by the scripts in `figures/`. They are run in order (Part01 to Part09) from the same working directory as the scripts in `preprocess/`, whose outputs make up most of their inputs. `TSS_manual_assignment.txt` (this folder) and `PlantTFDB_Ath_TF_binding_motifs.meme` (`data/external/`) should also be placed in the working directory, as well as the external files listed in the main README. Many of the outputs can be found in some form in the `data/` folder.

## Software used

- R v4.4.1
- bigWigAverageOverBed (version 2)
- streme, centrimo, sea and tomtom (MEME suite version 5.4.1)
- minimotif (version 1.3; [yamtk](https://github.com/bjmt/yamtk) commit b291730, including `scripts/to_bed.sh`)
- HOMER v4.11.1
