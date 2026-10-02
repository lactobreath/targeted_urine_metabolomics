# Targeted urine metabolomics

This repository contains the code for the targeted urine metabolomics data analysis published [here](toadd).

The repository contains a series of scripts to run the different analyses and produce the results presented in the manuscript. The analyses of the different scripts are briefly described below:

- `targeted_urine_metabolomics_with_sex.Rmd`: Processing of GC-MS data, differential analysis of urinary lactose metabolites by study group/genotype, longitudinal analysis by study group/genotype, intervention substrate and timepoints.

- `Comparison_by_covar.R`: exploration of potential covariate effects of urinary lactose metabolites levels.
- `Comparison_by_fluid_intake.R`: comparison of whole urine volume and fluid intake by study groups.
- `Multivariate_analysis_study_group_with_sex.R`: LDA of urinary lactose metabolite levels to discriminate study groups.
- `ROC_curves_metabolites_h2_3_pools_with_sex.R`: ROC curve analysis of urinary lactose metabolite levels to discriminate study groups.
- `Symptoms_correlation_final.R`: correlation analysis between urinary lactose metabolites and gastrointestinal symptoms.
