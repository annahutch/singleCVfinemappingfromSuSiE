# singleCVfinemappingfromSuSiE
Code to perform single CV fine-mapping with corrected coverage adjustment for regions with evidence of a single CV from SuSiE.

## 01_findregionstomallerfinemap.R

(NOTE NEED TO AMEND THIS SCRIPT FOR HOW MANY SUSIE SETS THERE ARE (i.e. cols set_1,...,set_3)

Usage: Rscript 01_findregionstomallerfinemap.R $file

where $file is the SuSiE output table from Chris

Output: regions2corrcov_*.RDS

where * is first 3 characters in $file

This is a list of data.frames of GWAS results for the regions where the Maller fine-mapping approach
with corrcoverage adjustment will be run (those regions with either (i) evidence of a single CV from SuSiE or
(ii) significant p-values but missing imputed values so SuSiE could not be run).

## 02_singleCVfinemapping.R

Usage: Rscript 02_singleCVfinemapping.R $file
within a slurm script as an array job

where $file is the SuSiE output table from Chris
and the array job is ran over how many regions to finemap (i.e. the length of the list from 01_findregionstomallerfinemap.R

Output: res_$file/ directory containing .RDS files with the output from single CV fine-mapping with corrcoverage adjustment

##03_collateallresults.R

Usage: Rscript 03_collateallresults.R $file

where $file is the SuSiE output table from Chris

Output: finalFMres_$file.RDS with the results from SuSiE and adjusted coverage fine-mapping

This final output file contains 2 lists. The first list contains sublists for the regions with evidence of multiple CVs and the results from fine-mapping with SuSiE. The second list contains sublists for the regions with evidence of a single CV and the results from single causal variant fine-mapping (Maller et al. 2012) with adjustment method (Hutchinson et al. 2020).



