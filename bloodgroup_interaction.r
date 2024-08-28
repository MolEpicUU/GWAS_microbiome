# load libraries
library(rio)
library(BiocParallel)
library(lmerTest)
library(metafor)

# linear and logistic mixed model function
mixed.fun <- function(formula, data, binomial = FALSE) {
  tryCatch(
    {
      warnings <- c()
      if (binomial) {
        withCallingHandlers(coef <- summary(glmer(formula, data = data, family = "binomial", nAGQ = 0))$coefficients[, c(1, 2, 4)], warning = function(x) warnings <<- c(warnings, x$message))
      } else {
        withCallingHandlers(coef <- summary(lmer(formula, data = data))$coefficients[, c(1, 2, 5)], warning = function(x) warnings <<- c(warnings, x$message))
      }
      coef <- as.data.frame(coef)
      colnames(coef) <- c("estimate", "se", "pval")
      coef$message <- paste(warnings, collapse = "; ")
      coef
    },
    error = function(e) {
      e$message
    }
  )
}

# run bloodgroup antigen interaction analysis
run.mixed <- function() {
  res <- lapply(1:nrow(loci), function(i) {
    locus <- loci$Locus[i]
    trait <- loci$Trait[i[
      model <- loci$Model[i]
      pheno <- pheno_bt
      if (model == "Linear") pheno <- pheno_qt
      binomial <- model == "Logistic"
      z <- cov[, c("age", "age_squared", "sex", "plate", paste0("PC", 1:10))]
      data <- data.frame(y = pheno[, trait], x = cov$abo_A, i = cov$secretor, r = cov$family, z)
      data <- data[which(complete.cases(data)), ]
      resa <- mixed.fun(as.formula(paste0("y ~ x * i + (1 | r) + ", paste(colnames(z), collapse = "+"))), data, binomial)
      data <- data.frame(y = pheno[, trait], x = cov$abo_B, i = cov$secretor, r = cov$family, z)
      data <- data[which(complete.cases(data)), ]
      resb <- mixed.fun(as.formula(paste0("y ~ x * i + (1 | r) + ", paste(colnames(z), collapse = "+"))), data, binomial)
      data <- data.frame(y = pheno[, trait], x = cov$lewis, i = cov$secretor, r = cov$family, z)
      data <- data[which(complete.cases(data)), ]
      resl <- mixed.fun(as.formula(paste0("y ~ x * i + (1 | r) + ", paste(colnames(z), collapse = "+"))), data, binomial)
      res <- as.data.frame(rbind(resa[nrow(resa), ], resb[nrow(resb), ], resl[nrow(resl), ]))
      res$trait <- trait
      res$locus <- locus
      res$antigen <- c("ABO A", "ABO B", "Lewis")
      res
  })
  do.call(rbind, res)
}

# import data and run interaction analysis
anno <- rio::import("scapis_metagenomics_mgs_annotations_v1.0.tsv")
loci <- rio::import("loci.tsv")

## scapis
cov <- rio::import("scapis_bloodgroup.tsv")
pheno_qt <- rio::import("phenotypes_quantitative_20231020.tsv")
pheno_bt <- rio::import("phenotypes_binary_20231012.tsv")
scapis <- run.mixed()

## mos
cov <- rio::import("mos_bloodgroup.tsv")
pheno_qt <- rio::import("phenotypes_quantitative_20240214.tsv")
pheno_bt <- rio::import("phenotypes_binary_20240214.tsv")
mos <- run.mixed()

## cosmc
cov <- rio::import("cosmc_bloodgroup.tsv")
pheno_qt <- rio::import("cosmc_phenotypes_quantitative_20231220.tsv")
pheno_bt <- rio::import("cosmc_phenotypes_binary_20231220.tsv")
cosmc <- run.mixed()

## smcc
cov <- rio::import("smcc_bloodgroup.tsv")
pheno_qt <- rio::import("smcc_phenotypes_quantitative_20231220.tsv")
pheno_bt <- rio::import("smcc_phenotypes_binary_20231220.tsv")
smcc <- run.mixed()

# perform meta-analysis
res <- lapply(1:nrow(scapis), function(i) {
  estimates <- c(scapis$estimate[i], mos$estimate[i], cosmc$estimate[i], smcc$estimate[i])
  standarderrors <- c(scapis$se[i], mos$se[i], cosmc$se[i], smcc$se[i])
  meta <- unlist(rma(yi = estimates, sei = standarderrors, method = "FE")[c("beta", "se", "pval", "I2")])
  meta <- as.data.frame(t(meta))
  colnames(meta) <- c("estimate", "se", "pval", "isq")
  meta$trait <- scapis$trait[i]
  meta$locus <- scapis$locus[i]
  meta$antigen <- scapis$antigen[i]
  meta
})
res <- do.call(rbind, res)
res$model <- loci$Model[match(res$trait, loci$Trait)]
res <- res[, c("locus", "trait", "antigen", "model", colnames(res)[which(!colnames(res) %in% c("locus", "trait", "antigen", "model"))])]
resa <- res[which(res$antigen == "ABO A"), which(!colnames(res) == "antigen")]
resb <- res[which(res$antigen == "ABO B"), c("estimate", "se", "pval")]
resl <- res[which(res$antigen == "Lewis"), c("estimate", "se", "pval")]
res <- data.frame(resa[, c("locus", "trait", "model", "estimate", "se", "pval")], resb, resl)

# export results
rio::export(res, "results/bloodgroup_interaction.tsv")

