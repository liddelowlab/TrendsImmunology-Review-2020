# Liddelow et al., 2020, REVIEW


#### Code used to re-analyze published single cell RNASeq (scRNAseq) astrocyte data included in:
Liddelow SA, Marsh SE, Stevens B. [Microglia and astrocytes in health and disease: Dynamic Duo or Partners in Crime?](citation)


## Code
Included is the code necessary to replicate the Seurat objects used for analysis.
The following versions of individual software/packages were used for analysis:
- R 3.6.2
- ComplexHeatmap 2.2.0
- cowplot 1.0.0
- data.table 1.12.8
- dplyr 0.8.5
- ggplot2 3.3.0
- psych 1.9.12.31
- purrr 0.3.4
- RColorBrewer 1.1-2
- scater 1.14.6
- scds 1.2.0
- scran 1.14.6
- Seurat 3.1.5
- SingleCellExperiment 1.8.0
- viridis 0.5.1 


Code included in repository:
- 200629_demultiplexing_mathys_TD.R - demultiplexing code prepared by Taitea Dukstra (TD)
- 200629_mathys_reanalysis_JSS.R - analysis code prepared by Jessica Sadick (JSS)


## Data
### Original Data
The data analyzed in this project can be directly downloaded from the following source:

- [Mathys et al., 2019 (Nature)](https://www.nature.com/articles/s41586-019-1195-2), PMID: [31042697](https://pubmed.ncbi.nlm.nih.gov/31042697/), Human, snRNAseq (10X 3' V2), [syn18485175](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn18485175)


## Acknowledgements
Data analysis on astrocytes was completed by Jessica Sadick and Taitea Dykstra.
Additional data analysis in this Review Article on myeloid lineage cells completed by Samueal Marsh (@samuel-marsh), and can be found at [Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts](https://github.com/samuel-marsh/Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts).

The analysis and results published here are in whole or in part are based on data obtained from the [AMP-AD Knowledge Portal](https://adknowledgeportal.synapse.org/). Samples for this study were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161, R01AG15819, R01AG17917, R01AG30146, R01AG36836, U01AG32984, U01AG46152, the Illinois Department of Public Health, and the Translational Genomics Research Institute.
Specific datasets for which additional analyses were performed are available from the Synapse database under accession number [syn18485175](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn18485175) (Mathys).  Additional ROSMAP data can be requested at [https://www.radc.rush.edu](https://www.radc.rush.edu)
 
