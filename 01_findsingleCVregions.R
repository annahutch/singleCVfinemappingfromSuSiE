# script to extract genomic regions with evidence of a single causal variant from SuSiE output

# input: 
# [1] file name of SuSiE output

# Rscript singleCV_finemapping.R file

library(dplyr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

####################################################################
####################################################################

# extract regions with evidence of a 
# single causal variant and write these to a file

file <- args # file with susie results
#blocks <- fread("blocks.txt")

res <- fread(args)

# limit data-frame to only regions that were fine-mapped by Chris (p.imp < 1e-6)
res_fm <- res[which(res$finemapped==T), ]

# split up into LD blocks
y <- split(res_fm, res_fm$block)

# find regions with at most one credible set
one_cs <- which(lapply(y, function(x) all(is.na(x[,c("set_2","set_3","set_4","set_5","set_6","set_7","set_8")]))) %>% unlist())

# exclude regions where there is no credible set
no_sets <- which(lapply(y[one_cs], function(x) all(is.na(x[,c("set_1")]))) %>% unlist())

if(length(no_sets) > 0) one_cs <- one_cs[-no_sets,]

saveRDS(y[one_cs], paste0("singleCVregions_", file, ".RDS"))

####################################################################
####################################################################