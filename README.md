# Small patches hold most rare species in fragmented landscapes

This repository contains the prepared community datasets and MATLAB code used to produce the figures and tables for the manuscript:

> **Small patches hold most rare species in fragmented landscapes**

The analysis examines whether rare species occur in small habitat patches because small patches are disproportionately attractive to rare species, or because small patches are numerically dominant in fragmented landscapes. It uses 15 community-abundance datasets from fragmented landscapes on five continents.

The repository includes prepared fragment-by-species abundance matrices, source from the LandFrag, FragSAD and Bovendorp datasets. These original datafiles are not reproduced

To create the key figures and tables from the manuscript, run the matlab scripts in the manuscript_figures/ folder.

## Directory structure

```text
.
├── README.md
├── data/
│   ├── study_metadata.csv
│   └── <study_id>.csv
│
└── manuscript_figures/
    ├── fig_main_synoptic.m
    ├── fig_main_perdataset.m
    ├── fig_supp_perdataset.m
    ├── run_setcover_and_table.m
    ├── map_study_locations.m
    └── output/
