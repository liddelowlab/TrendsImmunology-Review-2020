# Trends-Immunology-Review-2020


#### Code used to re-analyze published single cell RNASeq (scRNAseq) astrocyte data included in:
Liddelow SA, Marsh SE, Stevens B. [Microglia and astrocytes in health and disease: Dynamic Duo or Partners in Crime?](citation)


## Code
Included is the code necessary to replicate the Seurat objects used for analysis.
- Each R file specifies the version of Seurat used for analysis.
  - where multiple version of Seurat are used this has been noted.
  - Seurat V2.3.4 source package can be downloaded here from [CRAN Archive](https://cran.r-project.org/src/contrib/Archive/Seurat/) and installed from local source.

Code included in repository.
- 200629_demultiplexing_mathys_TD.R - demultiplexing code prepared by Taitea Dukstra (TD)
- 200629_mathys_reanalysis_JSS.R - analysis code prepared by Jessica Sadick (JSS)


## Data
### Original Data
The data analyzed in this project can be directly downloaded from the following source.

- [Mathys et al., 2019 (Nature)](https://www.nature.com/articles/s41586-019-1195-2), PMID: [31042697](https://pubmed.ncbi.nlm.nih.gov/31042697/), Human, snRNAseq (10X 3' V2), [syn18485175](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn18485175)


## Acknowledgements
Data analysis on astrocytes was completed by Jessica Sadick and Taite Dykstra .
Additional data analysis in this Review Article on myeloid lineage cells completed by Samueal Marsh (@samuel-marsh), and can be found at [Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts](https://github.com/samuel-marsh/Marsh_et-al_2020_scRNAseq_Dissociation_Artifacts).

The analysis and results published here are in whole or in part are based on data obtained from the [AMP-AD Knowledge Portal](https://adknowledgeportal.synapse.org/). Samples for this study were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161, R01AG15819, R01AG17917, R01AG30146, R01AG36836, U01AG32984, U01AG46152, the Illinois Department of Public Health, and the Translational Genomics Research Institute.
Specific datasets for which additional analyses were performed are available from the Synapse database under accession number [syn18485175](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage?Study=syn18485175) (Mathys).  Additional ROSMAP data can be requested at [https://www.radc.rush.edu](https://www.radc.rush.edu)
 
