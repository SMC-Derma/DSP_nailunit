# DSP_nailunit

+ Our previously conducted single-cell RNA analysis in human nail units confirmed the presence of a nail-specific mesenchymal cell population called onychofibroblasts within the onychodermis (Communications Biology volume 4, Article number: 692, 2021)

+ Here, we used a spatial transcriptomics technology and defined a cellular composition and spatially resolved expression profile of the human HF and nail unit. 
We also integrated spatially resolved expression data and single-cell RNA sequencing to investigate transcriptional similarity between two major skin appendages. 
To support this hypothesis, we estimated the absolute abundance of cell types of the nail unit using a SpatialDecon algorithm

+ This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq and DSP data for our preprint:
  + scRNA-seq of polydactyly can be accessed from the NCBI Gene Expression Omnibus database (accession code GSE158970). 
   
+ Version information:  
  + FASTQ reads were mapped to GRCh38.
  + Single-cell RNA-seq: Seurat package version 3.1.1 in R version 4.0.3 software.
  + Ligand-receptor interaction analysis was done with NicheNet ver.1.0
  + Cell abundance deconvolution was done with SpatialDecon ver. 1.0.0 
  + DEG analysis was done with edgeR ver 3.32.1      
