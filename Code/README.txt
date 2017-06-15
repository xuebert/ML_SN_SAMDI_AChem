This folder contains the code to generate results for Analytical Chemistry manuscript "Machine learning on SAMDI mass spectrometry signal to noise ratio improves peptide array designs"
Authors: Xue, Albert; Szymczak, Lindsey; Mrksich, Milan; Bagheri, Neda

To run any script, ensure that the current working directory is set to this Code folder

Figures 2, 4, 5, 6 and Supporting Information Figures 1, 2, 3, 4, and Supporting Tables 1 and 2 can all be created using scripts here
Significant post-processing was done to create final manuscript versions

######## description of subfolders #######
figure_dump contains all the pdf versions of figures (or tables)
figure_scripts converts the data or machine learning results into figures
ML_functions help with formatting different pieces of the machine learning
support_functions contain utility functions like color, loading data, and calculating p-values

The two machine learning scripts, predict_SN_script.R and sampling_peptides_script.R can be run as is to generate results.  Their corresponding figure scripts must be run after the machine learning scripts.

sampling_peptides_script.R was run on Northwestern HPC, Quest.  As a result, the sampling_peptides_results.RData file is in a different format than what would be created from this script.  The current make_figure6.R script works for the provided sampling_peptides_results.RData file but not for the output of the current sampling_peptides_script.R

Also, sampling_peptides_script will take a long time to run with many bootstraps.  It is recommended to parallelize the script.
