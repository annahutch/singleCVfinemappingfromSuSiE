# script to perform fine-mapping (with corrected coverage adjustment)
# on the genomic regions with evidence of a single CV or where 
# SuSiE cannot be run due to missing imputed values

# input: 
# [1] file name of SuSiE output

# submit as an array job with indices
# from 1 to the number of regions to finemap
# i.e. 1-(length of regions2corrcov_* file)

# need to be in parent directory of LD/

library(dplyr)
library(data.table)
library(corrcoverage)

args <- commandArgs(trailingOnly=TRUE)

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
i <- as.numeric(task_id_string)

####################################################################
####################################################################

file <- args[[1]] # file with susie results

y_fm <- readRDS(paste0("regions2corrcov_", substr(file, 1, 3), ".RDS"))

load(paste0("../LD/",names(y_fm)[i],".RData"))

# check 
if(length(MAF)!=dim(LD)[1]) print("MAF and LD diff sizes")

# match SNPs in GWAS to LD
x_subset <- y_fm[[i]][na.omit(match(names(MAF), y_fm[[i]]$snp)),]

# check 
if(length(MAF)!=dim(x_subset)[1]) {
  missing <- which(is.na(match(names(MAF), y_fm[[i]]$snp)))
  # remove this from MAF and LD objects
  MAF <- MAF[-missing]
  LD <- LD[-missing, -missing]
  print("MAF and x_subset diff sizes")
}

if(identical(x_subset$snp, names(MAF))==FALSE) print("Wrong order")

# calculate PPs
W = 0.2
z = x_subset$beta.orig/x_subset$se.orig
r = W^2/(W^2 + x_subset$se.orig^2)
pp = ppfunc(z, V = x_subset$se.orig^2)

# need to name these variables to list the variants in the credset in the output
mybhat <- x_subset$beta.orig
names(mybhat) <- x_subset$snp
names(pp) <- x_subset$snp

out <- tryCatch(
  {
    corrected_cs_bhat(bhat = mybhat,
                      V = x_subset$se.orig^2,
                      N0 = 0,  # N0 and N1 not actually usedd
                      N1 = 0,
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

#dir.create(paste0("res_", file))
saveRDS(list(out, x_subset), paste0("res/","res_",names(y_fm)[i],".RDS"))
