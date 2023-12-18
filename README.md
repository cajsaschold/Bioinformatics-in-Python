# Genetic Variations Analysis for Any Gene

## Project Overview

This project focuses on analyzing genetic variations within a specified genomic region for any gene. The data should be presented in Variant Call Format (VCF), and the analysis includes examining alternative allele frequencies, identifying nonsynonymous variants and their apperance in the specific genomic region. 

## Project Description

The data is from the 1000 Genomes project and is presented in the Variant Call Format (VCF), with genomes compared to a human reference genome GH19. The analysis is designed to be versatile and can be applied to any gene of interest. In this project, the script is used to analyze the PTPN22 gene within the genomic region 1:114,356,433-114,414,375 of chromosome 1. 

## Methods

The project involves Python version 3 for data analysis. Key methods include counting samples and variants, identifying nonsynonymous variants and their apperances in different populations, and creating histograms for alternative allele frequencies. 
The analysis is designed to be flexible, allowing users to specify the gene of interest.

## Code Execution

To run the analysis for any gene, execute the provided Python script:

python 1000genomes.py [vcf_file] [sample_info_file]
