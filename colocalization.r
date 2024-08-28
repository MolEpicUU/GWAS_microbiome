# load libraries
library(TwoSampleMR)
library(coloc)
library(rio)
library(BiocParallel)

# function to parse UKBB GWAS statistics
parse_ukbb <- function(ukbb, phenotype, type) {
  split <- strsplit(ukbb$variant, split = ":")
  split <- do.call(rbind, split)
  ukbb$chr <- split[, 1]
  ukbb$pos <- as.numeric(split[, 2])
  ukbb$oa <- split[, 3]
  ukbb$oa[which(ukbb$oa == ukbb$minor_allele)] <- split[, 4][which(ukbb$oa == ukbb$minor_allele)]
  ukbb$variant <- paste(ukbb$chr, ukbb$pos, sep = ":")
  ukbb$n <- median(ukbb$n_complete_samples, na.rm = TRUE)
  ukbb$type <- type
  ukbb$phenotype <- phenotype
  ukbb$id.phenotype <- phenotype
  ukbb <- ukbb[, c("type", "chr", "pos", "n", "variant", "beta", "se", "minor_allele", "oa", "minor_AF", "pval", "phenotype", "id.phenotype")]
  colnames(ukbb) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
  ukbb
}

# function to parse METAL GWAS statistics
parse_metal <- function(metal, phenotype, type) {
  split <- strsplit(metal$MarkerName, split = ":")
  split <- do.call(rbind, split)
  metal$chr <- split[, 1]
  metal$pos <- as.numeric(split[, 2])
  metal$variant <- paste(metal$chr, metal$pos, sep = ":")
  metal$n <- 16040
  metal$type <- type
  metal$phenotype <- phenotype
  metal$id.phenotype <- phenotype
  metal <- metal[, c("type", "chr", "pos", "n", "variant", "Effect", "StdErr", "Allele1", "Allele2", "Freq1", "P-value", "phenotype", "id.phenotype")]
  colnames(metal) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
  metal
}

# function to parse REGENIE GWAS statistics
parse_regenie <- function(regenie, phenotype, type) {
  regenie$variant <- paste(regenie$CHROM, regenie$GENPOS, sep = ":")
  regenie$n <- median(regenie$N, na.rm = TRUE)
  regenie$p <- 10^-regenie$LOG10P
  regenie$type <- type
  regenie$phenotype <- phenotype
  regenie$id.phenotype <- phenotype
  regenie <- regenie[, c("type", "CHROM", "GENPOS", "n", "variant", "BETA", "SE", "ALLELE1", "ALLELE0", "A1FREQ", "p", "phenotype", "id.phenotype")]
  colnames(regenie) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
  regenie
}

# colocalization function
coloc <- function(locus, gwas1, gwas2, chr, pos) {
  gwas1 <- gwas1[which(gwas1$chr == chr & abs(gwas1$pos - pos) < 500000), ]
  gwas1 <- gwas1[order(gwas1$pval), ]
  gwas1 <- gwas1[which(!duplicated(gwas1$SNP)), ]
  region1 <- gwas1[, -1:-4]
  region2 <- gwas2[, -1:-4]
  colnames(region1) <- c("SNP", "beta.exposure", "se.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "pval.exposure", "exposure", "id.exposure")
  colnames(region2) <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "pval.outcome", "outcome", "id.outcome")
  harmonized <- harmonise_data(region1, region2)
  harmonized <- harmonized[which(harmonized$mr_keep), ]
  harmonized <- harmonized[order(harmonized$pval.exposure), ]
  harmonized <- harmonized[which(!duplicated(harmonized$SNP)), ]
  split <- strsplit(harmonized$SNP, split = ":")
  split <- do.call(rbind, split)
  harmonized$pos <- as.numeric(split[, 2])
  region1 <- list(beta = harmonized$beta.exposure, varbeta = harmonized$se.exposure^2, snp = harmonized$SNP, position = harmonized$pos, type = gwas1$type[1], N = gwas1$n[1], MAF = harmonized$eaf.exposure)
  region2 <- list(beta = harmonized$beta.outcome, varbeta = harmonized$se.outcome^2, snp = harmonized$SNP, position = harmonized$pos, type = gwas2$type[1], N = gwas2$n[1], MAF = harmonized$eaf.outcome)
  pheno1 <- gwas1$phenotype[1]
  pheno2 <- gwas2$phenotype[1]
  p1 <- harmonized$pval.exposure
  p2 <- harmonized$pval.outcome
  warnings <- c()
  withCallingHandlers(out <- capture.output(res <- coloc.abf(dataset1 = region1, dataset2 = region2)), warning = function(x) warnings <<- c(warnings, x$message))
  save(locus, pheno1, pheno2, region1, region2, p1, p2, res, warnings, out, file = paste0("results/coloc_", tolower(gsub(",", "", locus)), "_", tolower(pheno1), "_", tolower(pheno2), ".rda"))
}

# clean glucose GWAS
gwas <- rio::import("raw/30740_irnt.gwas.imputed_v3.both_sexes.varorder.tsv")
gwas <- parse_ukbb(gwas, "Glucose", "quant")
rio::export(gwas, "clean/glucose.tsv.gz")

# clean IBD GWAS
gwas <- rio::import("raw/EUR.IBD.gwas_info03_filtered.assoc", format = "tsv")
gwas$variant <- paste(gwas$CHR, gwas$BP, sep = ":")
gwas$beta <- log(gwas$OR)
gwas$freq <- (gwas$FRQ_A_12882 * 12882 + gwas$FRQ_U_21770 * 21770) / (12882 + 21770)
gwas$n <- 86640
gwas$type <- "cc"
gwas$phenotype <- "IBD"
gwas$id.phenotype <- "IBD"
gwas <- gwas[, c("type", "CHR", "BP", "n", "variant", "beta", "SE", "A1", "A2", "freq", "P", "phenotype", "id.phenotype")]
colnames(gwas) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
rio::export(gwas, "clean/ibd.tsv.gz")

# clean stool frequency GWAS
gwas <- rio::import("raw/GCST90002250_buildGRCh37.tsv.gz")
gwas$variant <- paste(gwas$chromosome, gwas$base_pair_location, sep = ":")
gwas$n <- 167875
gwas$type <- "quant"
gwas$phenotype <- "Stool frequency"
gwas$id.phenotype <- "Stool frequency"
gwas <- gwas[, c("type", "chromosome", "base_pair_location", "n", "variant", "beta", "standard_error", "effect_allele", "other_allele", "effect_allele_frequency", "p_value", "phenotype", "id.phenotype")]
colnames(gwas) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
rio::export(gwas, "clean/stool_freq.tsv.gz")

# clean WHRadjBMI GWAS
gwas <- rio::import("raw/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz")
gwas$variant <- paste(gwas$CHR, gwas$POS, sep = ":")
gwas$n <- median(gwas$N, na.rm = TRUE)
gwas$type <- "quant"
gwas$phenotype <- "WHRadjBMI"
gwas$id.phenotype <- "WHRadjBMI"
gwas <- gwas[, c("type", "CHR", "POS", "n", "variant", "BETA", "SE", "Tested_Allele", "Other_Allele", "Freq_Tested_Allele", "P", "phenotype", "id.phenotype")]
colnames(gwas) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
rio::export(gwas, "clean/whradjbmi.tsv.gz")

# clean LDLC GWAS
gwas <- rio::import("raw/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz", format = "tsv")
gwas$pvalue <- as.numeric(gwas$pvalue)
gwas$variant <- paste(gwas$CHROM, gwas$POS_b37, sep = ":")
gwas$n <- median(gwas$N, na.rm = TRUE)
gwas$type <- "quant"
gwas$phenotype <- "LDLC"
gwas$id.phenotype <- "LDLC"
gwas <- gwas[, c("type", "CHROM", "POS_b37", "n", "variant", "EFFECT_SIZE", "SE", "ALT", "REF", "POOLED_ALT_AF", "pvalue", "phenotype", "id.phenotype")]
colnames(gwas) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
rio::export(gwas, "clean/ldlc.tsv.gz")

# clean SHBG GWAS
gwas <- rio::import("raw/GCST90012111_buildGRCh37.tsv.gz")
gwas$variant <- paste(gwas$chromosome, gwas$base_pair_location, sep = ":")
gwas$n <- 425097
gwas$type <- "quant"
gwas$phenotype <- "SHBG"
gwas$id.phenotype <- "SHBG"
gwas <- gwas[, c("type", "chromosome", "base_pair_location", "n", "variant", "beta", "standard_error", "effect_allele", "other_allele", "effect_allele_frequency", "p_value", "phenotype", "id.phenotype")]
colnames(gwas) <- c("type", "chr", "pos", "n", "SNP", "beta", "se", "effect_allele", "other_allele", "eaf", "pval", "phenotype", "id.phenotype")
rio::export(gwas, "clean/shbg.tsv.gz")

# import data
anno <- rio::import("scapis_metagenomics_mgs_annotations_v1.0.tsv")
loci <- rio::import("loci.tsv")

# perform colocalization for studywide microbial traits vs selected traits
BiocParallel::bplapply(1:nrow(loci), function(i) {
  locus <- loci$Locus[i]
  trait <- loci$Trait[i]
  chr <- loci$Chr[i]
  pos <- loci$Pos37[i]
  model <- loci$Model[i]
  file <- loci$File[i]
  type <- ifelse(model == "Linear", "quant", "cc")
  gwas1 <- rio::import(paste0("metal/output/", file))
  gwas1 <- parse_metal(gwas1, trait, type)
  files <- list.files("clean", full.names = TRUE)
  lapply(files, function(x) {
    gwas2 <- rio::import(x)
    coloc(locus, gwas1, gwas2, chr, pos)
  })
}, BPPARAM = BiocParallel::MulticoreParam(6))

# perform colocalization for studywide microbial traits vs studywide microbial traits 
BiocParallel::bplapply(1:nrow(loci), function(i) {
  locus <- loci$Locus[i]
  trait <- loci$Trait[i]
  chr <- loci$Chr[i]
  pos <- loci$Pos37[i]
  model <- loci$Model[i]
  file <- loci$File[i]
  type <- ifelse(model == "Linear", "quant", "cc")
  gwas1 <- rio::import(paste0("metal/output/", file))
  gwas1 <- parse_metal(gwas1, trait, type)
  lapply(loci$File[-i], function(x) {
    gwas2 <- rio::import(paste0("metal/output/", x))
    mgs_id <- gsub("qt_|bt_|alpha_|_1.txt", "", x)
    trait <- anno$maintax[match(mgs_id, anno$mgs_id)]
    if (mgs_id == "richness") trait <- "Richness"
    type <- ifelse(grepl("bt", x), "cc", "quant")
    gwas2 <- parse_metal(gwas2, trait, type)
    coloc(locus, gwas1, gwas2, chr, pos)
  })
}, BPPARAM = BiocParallel::MulticoreParam(6))

# perform colocalization for studywide microbial traits vs secondary bile acids
met_anno <- rio::import("scapis_metabolon_annotations_batchnorm_v2.1.tsv")
BiocParallel::bplapply(1:nrow(loci), function(i) {
  locus <- loci$Locus[i]
  trait <- loci$Trait[i]
  chr <- loci$Chr[i]
  pos <- loci$Pos37[i]
  model <- loci$Model[i]
  file <- loci$File[i]
  type <- ifelse(model == "Linear", "quant", "cc")
  gwas1 <- rio::import(paste0("metal/output/", file))
  gwas1 <- parse_metal(gwas1, trait, type)
  files <- list.files("sba/")
  id <- gsub("qt_|bt_|.gz", "", files)
  lapply(id, function(x) {
    file <- files[which(grepl(x, files))]
    gwas2 <- rio::import(paste0("sba/", file), format = "tsv")
    trait <- met_anno$CHEMICAL_NAME[match(x, met_anno$met_id)]
    type <- ifelse(grepl("bt", file), "cc", "quant")
    gwas2 <- parse_regenie(gwas2, trait, type)
    coloc(locus, gwas1, gwas2, chr, pos)
  })
}, BPPARAM = BiocParallel::MulticoreParam(6))

# make results table
files <- list.files("results", full.names = TRUE)
res <- lapply(files, function(x) {
  load(x)
  trait <- gsub("^.*/|.rda", "", x)
  trait <- unlist(strsplit(trait, split = "_"))
  res <- t(round(res$summary[2:6], 4))
  colnames(res) <- paste0("h", 0:4)
  data.frame(locus = locus, trait1 = pheno1, trait2 = pheno2, res)
})
res <- do.call(rbind, res)
res <- res[order(res$h4, decreasing = TRUE), ]
colnames(res) <- c("Locus", "Trait1", "Trait2", "H0", "H1", "H2", "H3", "H4")
res <- res[which(!res$Trait1 == res$Trait2), ]

# export results
rio::export(res, "results/coloc.tsv")

