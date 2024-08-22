## Spatial transcriptomics reveals structurally organized and distinct immune polarization in inflammatory cutaneous granulomatous disorders
GH repo written by Joe Daccache  
Updated on 04/02/2024  

The purpose of this GH repo is to provide the code from the 2024 publication using spatial transcriptomics on cutaneous granulomas disorders. The slices were obtained from the Yale biobank and processed using 10X Genomics' Visium assay. Libraries were sequenced by the Yale Genomic Core and analysis was performed using a combination of packages, namely Seurat and Giotto. If you have any questions, please feel free to reach out to myself at joseph.daccache@nyulangone.org  

### Analytical workflow
**Upstream Processing in spaceranger** 
*spaceranger*, the analysis tool developed by 10x Genomics was used to intiialy process our slides. The script can be found in `src/spaceranger.sh`

### wrapper functions
At the top of every processing script, the file containing wrapper functions `src/Spatial_wrappers.R` is read in. This file contains user defined functions used to perform the bulk of the analysis. The functions are genrally written to work on lists of seurat/giotto object, owing the nature of the study working on multiple objects at a time.  

### Initial processing in Seurat
These scripts were used to generate Figures 1B-E, 2B  
The preprocessing workflow, outlined in `src/preprocessing.R`, was performed on slices individually for qc purposes. Here the Seurat objects are created with filtering parameters and the user drags to remove spots not connected to rest of tissue space. Then the objects are clustered and annotated for histological region. Then the objects are filtered based on clustering/study criteria and more dragging any subsetted objects are then re-clustered and reannotated.  
  
Following qc and spot filtration, the user selects the optimal clustering resolution and annotates clusters. This entire process is repeated for the merging of samples within the same condition and is a very iterative process, like many single cell/spatial analyses end up being. The plots generated in this script are not in the final manuscript as they are mostly all preprocessing/quality control but they may be helpful to someone looking for specific code or ways to think about this data!  

### Initial processing in Giotto
Since we perform cell type deconvolution in Giotto, we need to convert Seurat objects to Giotto objets and construct spatial neighborhoods. Script is in `src/giotto_preprocessing.R`  

### Cell type deconvolution
Using *Giotto*, we performed cell type deconvolution of our visium spots/slides using previoulsy published scRNA-Seq datatsets generated in our lab. The script is memory heavy and should be run on a cluster. The script can be found in `src/celltype-deconvolution.R`  

### Module generation and enrichment  
For analyses presented in panel C of Figures 3-6, spatially correlated gene modules were generated across samples in the same condition using *Giotto*. The script for this analysis can be found in `src/module-generation-and-gsea.R`  

### Specific Visualizations
This script was used or modelled off of to generate Figures 1B&F, 2A-E, 3D-G, 4D-G, 5D-I, 6D-G, 7A&B  
For multiple figures, gene expression of specific genes was visualized in a variety of methods, the code to do so is included in `src/figures.R`, albeit some specific plots may be missing, the skeleton of the code is all there for any expression plot in the publication. Heatmaps of module expression *Figure panel A in figures 3-6* are generated in `srx/module-heatmaps.R`. The GSEA of macrophage marker comparisons are included (Figure 2E) in this script as well.  

### Spatial PCA analysis and plots
This script was used to generate figures in supplementals for each condition.  
This script, `src/spatial-pca.R` uses the *SpatialPCA* package to run spatially informed psuedotime trajectory inference, more or less ranking spots by their degree of similarity (in the umap space) to the granulomous spots. 
