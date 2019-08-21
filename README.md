# READ ME
This is the code repository for the paper *"ADHD symptoms map onto noise-driven structure-function decoupling between hub and peripheral brain regions"* available at Molecular Psychiatry here: [not yet available] and via preprint here: https://www.biorxiv.org/content/biorxiv/early/2019/04/11/606228.full.pdf

Authors: Luke J. Hearne*, Hsiang-Yuan Lin*, Paula Sanz-Leon, Wen-Yih Isaac Tseng, Susan Shur-Fen Gau, James A. Roberts, Luca Cocchi

![alt text](https://github.com/ljhearne/ADHDSCFC/blob/master/Figures/Figures-01.png "Fig. 1 Conceptual overview of the analysis pipeline. A. Analyses were conducted using a whole-brain parcellation including 214 cortical and subcortical regions. Replication analyses were performed using two alternative brain parcellations (see text). B. Structural (SC) and functional connectivity (FC) matrices were derived from diffusion spectrum imaging (DSI) and multi-echo resting-state fMRI data, respectively. Darker colors indicate higher normalized streamline counts (SC) and higher Fisher-z normalized Pearson’s correlation values between every possible pair of brain regions (FC). C. The topological organization of the SC matrices was examined to derive measures of different connection types: hub connections, feeder connections, and local connections. Individual-level correlations between SC and FC were used to estimate structure-function coupling, which was then analyzed with between-group statistics. D. A computational model was used to assess the potential neural mechanisms that lead to decreased structure-function coupling. Empirical SC was used as input in the model and model parameters were estimated by fitting to empirical FC. We systematically assessed if an increase in the noise heterogeneity in hub or peripheral nodes could result in a marked dissociation between functional and structural connectivity.")


## Contact
For information regarding the repository please contact luke.hearne[at]gmail.com otherwise, email the corresponding author from the papers.

## Instructions for coauthors
Folders should be organized like so:
```
project folder 
│
└───Data
│   │   AdultADHD_FS....xlsx
│   └───Schaefer114 (as Hsiang-Yuan organised)
│   │   │    ...
│   └───Schaefer214 (as Hsiang-Yuan organised)
│       │    ...
│   
└───Docs (this repository)
    │    ...
    └─── Results
        └─── Kx (Hub definition used, e.g., K15)
            └─── Atlas (e.g., 214)
                └─── Contains figures + a stats.txt file which contains statistical results.
    └─── Scripts
```
