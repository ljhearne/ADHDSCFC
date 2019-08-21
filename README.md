# READ ME
This is the code repository for the paper *"ADHD symptoms map onto noise-driven structure-function decoupling between hub and peripheral brain regions"* available at Molecular Psychiatry here: and via preprint here:

Authors: Luke J. Hearne, Hsiang-Yuan Lin, Paula Sanz-Leon, Wen-Yih Isaac Tseng, Susan Shur-Fen Gau, James A. Roberts, Luca Cocchi

## Abstract
Adults with childhood-onset attention-deficit hyperactivity disorder (ADHD) show altered whole-brain connectivity. However, the relationship between structural and functional brain abnormalities, the implications for the development of life-long debilitating symptoms, and the underlying mechanisms remain uncharted. We recruited a unique sample of 80 medication-naive adults with a clinical diagnosis of childhood-onset ADHD without psychiatric comorbidities, and 123 age-, sex-, and intelligence-matched healthy controls. Structural and functional connectivity matrices were derived from diffusion spectrum imaging and multi-echo resting-state functional MRI data. Hub, feeder, and local connections were defined using diffusion data. Individual-level measures of structural connectivity and structure-function coupling were used to contrast groups and link behavior to brain abnormalities. Computational modeling was used to test possible neural mechanisms underpinning observed group differences in the structure-function coupling. Structural connectivity did not significantly differ between groups but, relative to controls, ADHD showed a reduction in structure-function coupling in feeder connections linking hubs with peripheral regions. This abnormality involved connections linking fronto-parietal control systems with sensory networks. Crucially, lower structure-function coupling was associated with higher ADHD symptoms. Results from our computational model further suggest that the observed structure-function decoupling in ADHD is driven by heterogeneity in neural noise variability across brain regions. By highlighting a neural cause of a clinically meaningful breakdown in the structure-function relationship, our work provides novel information on the nature of chronic ADHD. The current results encourage future work assessing the genetic and neurobiological underpinnings of neural noise in ADHD, particularly in brain regions encompassed by fronto-parietal systems.
## Contact
luke.hearne[at]gmail.com
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
