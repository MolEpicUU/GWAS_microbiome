# load libraries
library(rio)

# create FUMA input
meta <- rio::import("results/meta.tsv")
split <- strsplit(meta$MarkerName, split = ":")
split <- do.call(rbind, split)
meta$marker <- paste(split[, 1], split[, 2], sep = ":")
meta$chr <- split[, 1]
meta$pos <- as.numeric(split[, 2])
meta$pval <- meta$"P-value"
meta <- meta[order(meta$pval), ]
meta$rs_id[which(!grepl("rs", meta$rs_id))] <- NA
meta <- meta[which(!is.na(meta$rs_id)), ]
meta <- meta[which(!duplicated(meta$marker)), ]
meta <- meta[which(!duplicated(meta$rs_id)), ]
meta <- meta[which(meta$chr != "23"), ]

# export
rio::export(meta[, c("rs_id", "chr", "pos", "pval")], "fuma.tsv")
