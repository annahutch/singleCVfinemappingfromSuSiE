# script to perform fine-mapping (with corrected coverage adjustment)
# on the genomic regions with evidence of a single CV

# input: 
# [1] file name of SuSiE output
# [2] number of cases in GWAS
# [3] number of controls in GWAS

# submit as an array job with indices
# from 1 to the number of regions with 
# evidence of a single signal
# i.e. 1-(length of singleCVregions_* file)

# need to be in parent directory of mafld/

library(dplyr)
library(data.table)
library(corrcoverage)

args <- commandArgs(trailingOnly=TRUE)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(task_id_string)

####################################################################
####################################################################

file <- args[1] # file with susie results

y_fm <- readRDS(paste0("singleCVregions_", file, ".RDS"))

load(paste0("../mafld/tmp_mafld_",names(y_fm)[i],".RData"))

# check 
if(length(MAF)!=dim(LD)[1]) print("MAF and LD diff sizes")

# match SNPs
x_subset <- x[na.omit(match(names(MAF), x$snp)),]

# check 
if(length(MAF)!=dim(x_subset)[1]) {
  missing <- which(is.na(match(names(MAF), x$snp)))
  # remove this from MAF and LD objects
  MAF <- MAF[-missing]
  LD <- LD[-missing, -missing]
  print("MAF and x_subset diff sizes")
}

if(identical(x_subset$snp, names(MAF))==FALSE) print("Wrong order")

N0 = args[3]
N1 = args[2]
N = N0 + N1

# calculate PPs
W = 0.2
z = x_subset$beta.imp/x_subset$se.imp
r = W^2/(W^2 + x_subset$se.imp^2)
pp = ppfunc(z, V = x_subset$se.imp^2)

# need to name these variables to list the variants in the credset in the output
mybhat <- x_subset$beta.imp
names(mybhat) <- x_subset$snp
names(pp) <- x_subset$snp

out <- tryCatch(
  {
    corrected_cs_bhat(bhat = mybhat,
                      V = x_subset$se.imp^2,
                      N0 = N0, 
                      N1 = N1,
                      Sigma = LD,
                      lower = 0.5,
                      desired.cov = 0.95,
                      acc = 0.001,
                      max.iter = 50)
  },
  error=function(cond) {
    new.out <- credset(pp, thr = 0.95)
    new.out$credset <- names(MAF)[new.out$credset]
    return(new.out)
  }
)

# add column to data.frame with my PPs and credset
x_subset$PP <- pp
x_subset$corr_credset <- FALSE
x_subset$corr_credset[match(out$credset, x_subset$snp)] <- TRUE

dir.create(paste0("res_", file))
saveRDS(list(out, x_subset), paste0("res_", file, "/res_",names(y_fm)[i],".RDS"))