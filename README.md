# NT_GSL_analysis

Repository for the publication "Unraveling the Glycosphingolipid Metabolism Leveraging Transcriptome-weighted Network Analysis on Neuroblastic Tumors" by 

Arsenij Ustjanzew (1*), Annekathrin Silvia Nedwed (1), Federico Marini (1, 2), Roger Sandhoff (3, 4), Jörg Faber (5, 6), and Claudia Paret (5, 6)

1 Institute of Medical Biostatistics, Epidemiology and Informatics (IMBEI), University Medical Center of the Johannes Gutenberg-University Mainz, 55131, Mainz, Germany

2 Research Center for Immunotherapy (FZI), Langenbeckstraße 1, 55131 Mainz, Germany

3 Lipid Pathobiochemistry, German Cancer Research Center, 69120 Heidelberg, Germany

4 Helmholtz-Institute for Translational Oncology Mainz (HI-TRON), 55131 Mainz, Germany

5 Department of Pediatric Hematology/Oncology, Center for Pediatric and Adolescent Medicine, University Medical Center of the Johannes Gutenberg-University Mainz, 55131, Mainz, Germany

6 University Cancer Center (UCT), University Medical Center of the Johannes Gutenberg-University Mainz, 55131 Mainz, Germany

*	Correspondence: arsenij.ustjanzew@uni-mainz.de;

# Overview

This repository contains following files:

* substrate_graph.Rds - contains the used GSL network as igraph object
* generating_network.Rmd - contains the steps of GSL network generation
* read_xenabrowser_data.R - Code for reading the TCGA-GTEx-TARGET data, creating a SummarizedExperiment object, and subsetting it to TARGET
* functions.R - contain functions for: 1) Assigning expression values to edges of an igraph object, 2) computing the reaction activity of a gene expression matrix, 3) computing transition probabilities with different adjustment options
* NB_GNB_pathanalysis.Rmd - contains the full analysis with visualizations of the TARGET NB/GNB dataset
* NB_GN_pathanalysis.Rmd - contains the full analysis with visualizations of the GSE147635 NB/GN dataset
* NB_GN_FC_image.Rmd - Code for creating figures of the Ganglioside graph and fold-change dotplot (Neuroblastoma vs. Ganglioneuroma)
* NB_GNB_FC_image.Rmd - Code for creating figures of the Ganglioside graph and fold-change dotplot (Neuroblastoma vs. Ganglioneuroblastoma)
* NBmycn_GNB_FC_image.Rmd - Code for creating figures of the Ganglioside graph and fold-change dotplot (Neuroblastoma MYCN amplified vs. Ganglioneuroblastoma)
* NBmycn_NB_FC_image.Rmd - Code for creating figures of the Ganglioside graph and fold-change dotplot (Neuroblastoma MYCN amplified vs. Neuroblastoma without MYCN amplification)

# Notes

This repository does not contain the used data (input data of the Rmd scripts) nor the output images/ files generated during the analysis.
* The TARGET NB/GNB dataset has to be downloaded and preprocessed from [here](https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) ("A combined cohort of TCGA, TARGET and GTEx samples"). The used file were the "gene expression RNAseq - RSEM expected_count" ("TcgaTargetGtex_gene_expected_count") and the "phenotype - TCGA TARGET GTEX selected phenotypes" ("TcgaTargetGTEX_phenotype.txt")
These data was combined to a SummarizedExperiment class object and subsetted to Neuroblastoma samples.
* The raw fastq files of the GSE147635 dataset were downloaded from Gene Expression Omnibus and processed as described in the manuscript.
* The R Markdown files are not intended to be knitted to a document and are meant to represent steps of the analysis.  


