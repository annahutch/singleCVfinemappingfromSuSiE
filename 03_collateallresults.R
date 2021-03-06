# script to collate the results from susie fine-mapping
# and adjusted credible set fine-mapping for single CV regions

# input: file name of SuSiE output

library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
file <- args # file with susie results

susie_res <- fread(args)

# limit data-frame to only regions that were fine-mapped by Chris (p.imp < 1e-6)
susie_res_fm <- susie_res[which(susie_res$finemapped==T), ]

# split up into LD blocks
susie_res_list <- split(susie_res_fm, susie_res_fm$block)

####################

singleCVfm_files <- list.files(path = "res", pattern=".RDS", full=TRUE)

singleCVfm_list <- lapply(singleCVfm_files, readRDS)

# extract second element from each list which is the full results
singleCVfm_shortlist <- sapply(singleCVfm_list, "[[", 2, simplify = F)

# remove from final list element 1 the single CV regions

all_regions <- lapply(susie_res_list, function(x) x$block[1]) %>% unlist()
singleCV_regions <- lapply(singleCVfm_shortlist, function(x) x$block[1]) %>% unlist()

susie_res_only <- susie_res_list[-na.omit(match(singleCV_regions, all_regions))]

final_fm_res <- list(susie_res_only, singleCVfm_shortlist)
                           
#length(final_fm_res[[1]])
#length(final_fm_res[[2]])
#x <- readRDS("regions2corrcov_UC_.RDS")
                           
saveRDS(final_fm_res, paste0("finalFMres_", file, ".RDS"))
