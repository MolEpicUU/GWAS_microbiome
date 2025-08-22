# centered log ratio transformation
clr <- function(data) {
    zeroes <- data == 0
    data <- data + min(data[data > 0])
    data <- log(data) - rowMeans(log(data))
    data[zeroes] <- -Inf
    data
}

# rank inverse normal transformation
rin <- function(x) qnorm((rank(x, "keep") - 0.5) / sum(!is.na(x)))

# generate relative abundances for higher taxa from species
mgs <- rio::import("scapis_metagenomics_mgs_relative_abundances_v1.0.tsv")
anno <- rio::import("scapis_metagenomics_mgs_annotations_v1.0.tsv")
merge_tax <- function(level, anno, mgs) {
    merged <- bplapply(unique(anno[, level]), function(x) {
        mgs <- mgs[, which(colnames(mgs) %in% anno$mgs_id[which(anno[, level] == x)]), drop = FALSE]
        rowSums(mgs)
    }, BPPARAM = MulticoreParam(16))
    merged <- do.call(cbind, merged)
    colnames(merged) <- paste0(level, "___", gsub(" ", "_", unique(anno[, level])))
    merged
}
anno <- anno[which(anno$mgs_id %in% colnames(mgs)), ]
taxa <- lapply(c("genus", "family", "order", "class", "phylum", "superkingdom"), function(x) merge_tax(x, anno, mgs))
taxa <- do.call(cbind, taxa)

# generate gm modules
mgs <- rio::import("scapis_metagenomics_mgs_relative_abundances_v1.0.tsv")
gm_anno <- rio::import("scapis_metagenomics_gm_annotations_v1.0.tsv")
gm <- lapply(unique(gm_anno$gm_id), function(x) {
    mgs_ids <- unlist(strsplit(gm_anno$mgs_id[which(gm_anno$gm_id == x)], split = ","))
    rowSums(mgs[, which(colnames(mgs) %in% mgs_ids)])
})
gm <- do.call(cbind, gm)
colnames(gm) <- unique(gm_anno$gm_id)
