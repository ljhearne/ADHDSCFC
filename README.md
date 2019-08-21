# READ ME
This is the code repository for the paper *"ADHD symptoms map onto noise-driven structure-function decoupling between hub and peripheral brain regions"* available at Molecular Psychiatry here: [not yet available] and via preprint here: https://www.biorxiv.org/content/biorxiv/early/2019/04/11/606228.full.pdf

Authors: Luke J. Hearne*, Hsiang-Yuan Lin*, Paula Sanz-Leon, Wen-Yih Isaac Tseng, Susan Shur-Fen Gau, James A. Roberts, Luca Cocchi    

<img src="https://github.com/ljhearne/ADHDSCFC/blob/master/Figures/Figures-01.png" width="500"/><br/>
**Fig. 1** Conceptual overview of the analysis pipeline.<br/><br/><br/>
<img src="https://github.com/ljhearne/ADHDSCFC/blob/master/Figures/Figures-02.png" width="500"/><br/>
**Fig. 2** Structure-function relationships in drug-naïve adults with ADHD and healthy matched controls<br/><br/><br/>
<img src="https://github.com/ljhearne/ADHDSCFC/blob/master/Figures/Figures-03.png" width="500"/><br/>
**Fig. 3** Modeling the effect of noise heteroscedasticity on structure-function coupling. <br/><br/><br/>

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

