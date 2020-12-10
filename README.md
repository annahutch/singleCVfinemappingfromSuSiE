# singleCVfinemappingfromSuSiE
Code to perform single CV fine-mapping with corrected coverage adjustment for regions with evidence of a single CV from SuSiE.

## 01_findsingleCVregions.R

Usage: Rscript 01_findsingleCVregions.R $file

where $file is the SuSiE output table from Chris

Output: finalFMres_$file.RDS

which contains a list of genomic regions with evidence of a single CV from SuSiE

## 02_singleCVfinemapping.R

Usage: Rscript 02_singleCVfinemapping.R $file $N1 $N0 
within a slurm script as an array job

where $file is the SuSiE output table from Chris, $N1 is the number of cases in the GWAS and $N0 is the number of controls in the GWAS
and the array job is ran over how many regions have evidence of a single CV 

Output: res_$file/ directory containing .RDS files with the output from single CV fine-mapping with corrcoverage adjustment for the regions with evidence of a single CV from SuSiE

##03_collateallresults.R

Usage: Rscript 03_collateallresults.R $file

where $file is the SuSiE output table from Chris

Output: finalFMres_$file.RDS with the results from SuSiE and adjusted coverage fine-mapping

This final output file contains 2 lists. The first list contains sublists for the regions with evidence of multiple CVs and the results from fine-mapping with SuSiE. The second list contains sublists for the regions with evidence of a single CV and the results from single causal variant fine-mapping (Maller et al. 2012) with adjustment method (Hutchinson et al. 2020).



