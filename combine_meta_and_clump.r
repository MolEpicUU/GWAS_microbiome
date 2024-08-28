# load libraries
library(rio)
library(BiocParallel)
library(meta)

# import data
gwas_anno <- rio::import("merged.tsv")

# generate meta-analysis results table
files <- list.files("metal/output", full.names = TRUE)
files <- files[which(!grepl("info", files))]
meta <- BiocParallel::bplapply(files, function(x) {
  data <- rio::import(x, format = "tsv")
  lambda <- median(qchisq(data$"P-value", 1, lower.tail = FALSE)) / qchisq(0.5, 1)
  data <- data[order(data$"P-value"), ]
  data2 <- data[1:10000, ]
  data3 <- data[10001:nrow(data), ]
  data3 <- data3[seq(1, nrow(data3), length.out = 20000), ]
  data <- rbind(data2, data3)
  x <- gsub("^.*/", "", x)
  x <- gsub("_1.txt", "", x)
  data$model <- ifelse(gsub("_.*$", "", x) %in% c("qt", "alpha"), "linear", "logistic")
  data$trait <- gsub("^.*_", "", x)
  data$lambda <- lambda
  data
}, BPPARAM = BiocParallel::MulticoreParam(16))
meta <- do.call(rbind, meta)
meta <- meta[which(!grepl("[?]", meta$Direction)), ]
meta$rs_id <- gwas_anno$ID[match(meta$MarkerName, gwas_anno$MarkerName)]

# export table for manhattan plot
rio::export(meta[, c("MarkerName", "rs_id", "P-value", "trait", "model", "lambda")], "results/meta_manhattan.tsv")

# genome-wide and study-wide results
meta <- meta[which(meta$"P-value" < 5e-8), ]

# add Isq CI
files <- list.files("metal/data", full.names = TRUE, pattern = , recursive = TRUE)
ci <- BiocParallel::bplapply(unique(meta$trait), function(x) {
  meta <- meta[which(meta$trait == x), ]
  data <- lapply(files[which(grepl(x, files))], function(x) {
    x <- rio::import(x, format = "tsv")
    x[match(meta$MarkerName, x$ID), c("BETA", "SE")]
  })
  data <- do.call(cbind, data)
  res <- apply(data, 1, function(x) {
    unlist(meta::metagen(x[c(1, 3, 5, 7)], x[c(2, 4, 6, 8)], random = FALSE)[c("I2", "lower.I2", "upper.I2")])
  })
  res <- as.data.frame(t(res))
  res$trait <- x
  res$MarkerName <- meta$MarkerName
  res
}, BPPARAM = BiocParallel::MulticoreParam(8))
ci <- do.call(rbind, ci)
meta$IsqLower <- ci$lower.I2[match(paste(meta$trait, meta$MarkerName), paste(ci$trait, ci$MarkerName)]
meta$IsqUpper <- ci$upper.I2[match(paste(meta$trait, meta$MarkerName), paste(ci$trait, ci$MarkerName)]

# export meta-analysis table
rio::export(meta, "results/meta.tsv")

# split in alpha diversity and species traits
alpha <- meta[which(meta$trait %in% c("richness", "shannon", "invsimpson")), ]
mgs <- meta[which(!meta$trait %in% c("richness", "shannon", "invsimpson")), ]

# clump genome-wide species
mgs <- mgs[which(mgs$"P-value" < (5e-8)),]
lapply(unique(mgs$trait), function(x) {
  mgs <- mgs[which(mgs$trait == x), ]
  split <- strsplit(mgs$MarkerName, split = ":")
  split <- do.call(rbind, split)
  mgs <- data.frame(SNP = mgs$rs_id, P = mgs$"P-value")
  chr <- split[, 1]
  chr[which(chr == "23")] <- "X"
  lapply(unique(chr), function(y) {
    mgs[which(chr == y), ]
    rio::export(mgs, paste0("clump/", x, "_", y, "_genomewide.tsv"))
    system(paste0("plink2 --pfile SCAPIS-GENETICS-PETITION-507-IMPUTED.EurQC_", y, " --clump clump/", x, "_", y, "_genomewide.tsv --out clump/", x, "_", y, "_genomewide --clump-r2 0.001 --clump-kb 10000"))
  })
})

# clump genome-wide alpha diversity traits
alpha <- alpha[which(alpha$"P-value" < (5e-8)), ]
lapply(unique(alpha$trait), function(x) {
  alpha <- alpha[which(alpha$trait == x), ]
  split <- strsplit(alpha$MarkerName, split = ":")
  split <- do.call(rbind, split)
  alpha <- data.frame(SNP = alpha$rs_id, P = alpha$"P-value")
  chr <- split[, 1]
  chr[which(chr == "23")] <- "X"
  lapply(unique(chr), function(y) {
    alpha[which(chr == y), ]
    rio::export(alpha, paste0("clump/", x, "_", y, "_genomewide.tsv"))
    system(paste0("plink2 --pfile SCAPIS-GENETICS-PETITION-507-IMPUTED.EurQC_", y, " --clump clump/", x, "_", y, "_genomewide.tsv --out clump/", x, "_", y, "_genomewide --clump-r2 0.001 --clump-kb 10000"))
  })
})

# clump study-wide species
mgs <- mgs[which(mgs$"P-value" < (5e-8 / 921)), ]
lapply(unique(mgs$trait), function(x) {
  mgs <- mgs[which(mgs$trait == x), ]
  split <- strsplit(mgs$MarkerName, split = ":")
  split <- do.call(rbind, split)
  mgs <- data.frame(SNP = mgs$rs_id, P = mgs$"P-value")
  chr <- split[, 1]
  chr[which(chr == "23")] <- "X"
  lapply(unique(chr), function(y) {
    mgs[which(chr == y), ]
    rio::export(mgs, paste0("clump/", x, "_", y, "_studywide.tsv"))
    system(paste0("plink2 --pfile SCAPIS-GENETICS-PETITION-507-IMPUTED.EurQC_", y, " --clump clump/", x, "_", y, "_studywide.tsv --out clump/", x, "_", y, "_studywide --clump-r2 0.001 --clump-kb 10000"))
  })
})

# clump study-wide alpha diversity traits
alpha <- alpha[which(alpha$"P-value" < (5e-8 / 3)), ]
lapply(unique(alpha$trait), function(x) {
  alpha <- alpha[which(alpha$trait == x), ]
  split <- strsplit(alpha$MarkerName, split = ":")
  split <- do.call(rbind, split)
  alpha <- data.frame(SNP = alpha$rs_id, P = alpha$"P-value")
  chr <- split[, 1]
  chr[which(chr == "23")] <- "X"
  lapply(unique(chr), function(y) {
    alpha[which(chr == y), ]
    rio::export(alpha, paste0("clump/", x, "_", y, "_studywide.tsv"))
    system(paste0("plink2 --pfile SCAPIS-GENETICS-PETITION-507-IMPUTED.EurQC_", y, " --clump clump/", x, "_", y, "_studywide.tsv --out clump/", x, "_", y, "_studywide --clump-r2 0.001 --clump-kb 10000"))
  })
})

# make results table
files <- list.files("clump/", full.names = TRUE, pattern = "clumps")
files <- files[which(!grepl("missing", files))]
clump <- lapply(files, function(x) {
  data <- rio::import(x, format = "tsv")
  split <- strsplit(gsub("^.*/|.clumps", "", x), split = "_")
  split <- do.call(rbind, split)
  data.frame(trait = split[, 1], level = split[, 3], data, check.names = FALSE)
})
clump <- do.call(rbind, clump)
rio::export(clump, "results/clump.tsv")
