# extract genomic regions with evidence of a single causal variant from SuSiE output and
# significant regions (p<5*10^-8) where SuSiE could not be run due to missing imputed values

# input: 
# [1] file name of SuSiE output

library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

####################################################################
####################################################################

# extract regions with evidence of a single causal variant

file <- args # file with susie results
#blocks <- fread("blocks.txt")

res <- fread(args)

# limit data-frame to only regions that were fine-mapped by Chris (p.imp < 1e-6)
res_fm <- res[which(res$finemapped==T), ]

# split up into LD blocks
y <- split(res_fm, res_fm$block)

# find regions with at most one credible set
one_cs <- which(lapply(y, function(x) all(is.na(x[,c("set_2","set_3")]))) %>% unlist())

# exclude regions where there is no credible set
no_sets <- which(lapply(y[one_cs], function(x) all(is.na(x[,c("set_1")]))) %>% unlist())

if(length(no_sets) > 0) one_cs <- one_cs[-no_sets]

####################################################################
####################################################################

# find significant regions where SuSiE could not be run
# due to missing imputed values

res_split <- split(res, res$block)

# find regions that do not have p.imp but do have significant p
sig <- lapply(res_split, function(x)  all(is.na(x$p.imp)) & min(x$p) < 5*10^-8 ) %>% unlist()

sig_res <- which(sig==TRUE)

regions_to_corrcov <- c(y[one_cs], res_split[sig_res])

saveRDS(regions_to_corrcov, paste0("regions2corrcov_", substr(file, 1, 3), ".RDS"))