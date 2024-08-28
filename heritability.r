# load libraries
library(rio)
library(BiocParallel)

# import data
anno <- rio::import("scapis_metagenomics_mgs_annotations_v1.0.tsv")
gwas_anno <- rio::import("merged.tsv")

# calculate GCTA tagging per chromosome
BiocParallel::bplapply(c(1:22, "X"), function(x) {
  sam <- rio::import(paste0("SCAPIS-GENETICS-PETITION-507-IMPUTED.EurQC_", x , ".psam"), format = "tsv")
  sam <- sample(sam$IID, 2000)
  sam <- data.frame(sam, sam)
  rio::export(sam, paste0(x, "_sam.txt"), col.names = FALSE)
  var <- rio::import(paste0("SCAPIS-GENETICS-PETITION-507-IMPUTED.EurQC_", x , ".pvar"), format = "tsv")
  info <- strsplit(var$INFO, split = ";")
  info <- sapply(info, function(x) grep("INFO=", x, value = TRUE))
  info <- as.numeric(gsub("INFO=", "", info))
  maf <- strsplit(var$INFO, split = ";")
  maf <- sapply(maf, function(x) grep("RefPanelAF=", x, value = TRUE))
  maf <- as.numeric(gsub("RefPanelAF=", "", info))
  var <- var[order(maf, decreasing = TRUE), ]
  var <- var[which(info > 0.8 & !var$ID %in% var$ID[which(duplicated(var$ID))] & !var$POS %in% var$POS[which(duplicated(var$POS))]), ]
  rio::export(as.data.frame(var$ID), paste0(x, "_var.txt"), col.names = FALSE)
  system(paste0("plink2 --pfile SCAPIS-GENETICS-PETITION-507-IMPUTED.EurQC_", x, " --make-bed --maf 0.01 --extract ", paste0(x, "_var.txt"), " --keep ", paste0(x, "_sam.txt"), " --out ", x))
  system(paste0("/ldak6.linux --calc-tagging gcta_", x, " --bfile ", x, " --power -1"))
  file.remove(paste0(x, ".bed"))
  file.remove(paste0(x, ".bim"))
  file.remove(paste0(x, ".fam"))
  file.remove(paste0(x, "_var.txt"))
  file.remove(paste0(x, "_sam.txt"))
}, BPPARAM = BiocParallel::MulticoreParam(6))

# merge tagging files
files <- paste0("gcta_", c(1:22, "X"), ".tagging")
rio::export(as.data.frame(files), "gcta.txt", col.names = FALSE)
system(paste0("ldak6.linux --join-tagging gcta --taglist gcta.txt"))

# estimate heritabilities
files <- list.files("metal/output/")
files <- files[which(!grepl("info", files))]
BiocParallel::bplapply(files, function(x) {
  data <- rio::import(paste0("metal/output/", x))
  data$A1 <- toupper(data$Allele1)
  data$A2 <- toupper(data$Allele2)
  data$Predictor <- gwas_anno$ID[match(data$MarkerName, gwas_anno$MarkerName)]
  data <- data[which(!is.na(data$Predictor)), ]
  data <- data[, c("Predictor", "A1", "A2", "Effect", "P-value")]
  colnames(data) <- c("Predictor", "A1", "A2", "Direction", "P")
  data <- data[order(data$P), ]
  data <- data[which(!duplicated(data$Predictor)), ]
  data$n <- 16040
  rio::export(data, paste0("clean_", x))
  trait <- gsub("alpha_|qt_|bt_|_1.txt", "", x)
  system(paste0("ldak6.linux --sum-hers gcta_", trait, " --summary ", paste0("clean_", x), " --tagfile gcta.tagging --check-sums NO --cutoff 0.01"))
  file.remove(paste0("clean_", x))
}, BPPARAM = BiocParallel::MulticoreParam(6))

# make results table
traits <- gsub("alpha_|qt_|bt_|_1.txt", "", files)
data <- lapply(traits, function(x) {
  gcta <- readLines(paste0("gcta_", x, ".hers"))
  gcta <- gcta[length(gcta)]
  gcta <- strsplit(gcta, split = " ")
  gcta <- unlist(gcta)
  heritability <- gcta[2]
  sd <- gcta[3]
  data.frame(trait = x, heritability, sd)
})
data <- do.call(rbind, data)
data$maintax <- anno$maintax[match(data$trait, anno$mgs_id)]
data$maintax[which(is.na(data$maintax))] <- sapply(data$trait[which(is.na(data$maintax))], function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))))
data <- data[, c("maintax", "heritability", "sd")]
colnames(data) <- c("Trait", "Heritability", "SD")
data <- data[order(data$Heritability, decreasing = TRUE), ]

# export results
rio::export(data, "results/heritability.tsv")
