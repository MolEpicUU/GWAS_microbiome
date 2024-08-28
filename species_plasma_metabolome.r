# load libraries
library(rio)
library(ppcor)

# partial Spearman model
spearman.fun <- function(y) {
  res <- BiocParallel::bplapply(colnames(met), function(x) {
    tryCatch({
      data <- data.frame(y, x = met[, x], pheno)
      data <- data[complete.cases(data), ]
      data <- data[, apply(data, 2, function(x) length(unique(x)) > 1)]
      data <- as.data.frame(model.matrix(~., data)[, -1])
      warnings <- c()
      withCallingHandlers(res <- ppcor::pcor.test(data$y, data$x, data[, which(!colnames(data) %in% c("y", "x"))], method = "spearman"), warning = function(x) warnings <<- c(warnings, x$message))
      data.frame(met_id = x, rho = res$estimate, pval = res$p.value, n = res$n, message = paste(warnings, collapse = "; "))
    }, error = function(e) {
      data.frame(met_id = x, rho = NA, pval = NA, n = nrow(data), message = paste(e$message, collapse = "; "))
    })
  }, BPPARAM = BiocParallel::MulticoreParam(16))
  do.call(rbind, res)
}

# import data
load("data.rda")
mgs_ids <- rio::import("mgs_ids.tsv", header = FALSE)[, 1]
pheno <- pheno[, which(!colnames(pheno) %in% c("scapis_id", "sample", "extraction_plate"))]

# run partial Spearman models
res <- lapply(mgs_ids, function(x) {
  res <- spearman.fun(abun[, mgs_id])
  data.frame(mgs_id = x, res)
})
res <- do.call(rbind, res)

# export results
rio::export(res, "results/plasma_metabolome.tsv")
