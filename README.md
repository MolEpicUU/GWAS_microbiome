# Analysis code for "Genome-wide association study highlights the role of the intestinal molecular environment in human microbiota variation" by Dekkers et al.

The GWAS pipeline used to generate the cohort-specific GWAS results can be found at: https://github.com/MolEpicUU/GWAS_scripts

# Scripts

* meta_analysis.r is the R-script used to meta-analyse the cohort-specific GWAS results
* heritability.r is the R-script used to calculate heritability estimates for the meta-analysis results
* combine_meta_and_clump.r is the R-script used to combine the meta-analysis results and perform LD clumping
* prepare_fuma.r is the R-script used to create input for FUMA analysis of the meta-analysis results
* colocalization.r is the R-script used for the colocalization analyses
* bloodgroup_interaction.r is the R-script used to perform antigen - secretor status interaction analysis
* species_plasma_metabolome.r is the R-script used to perform microbioal species plasma metabolite analysis
* stool_metabolome.r is the R-script used to lookup SNP - stool metabolite associations
