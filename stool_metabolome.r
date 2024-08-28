# load libraries
library(rio)
library(biomaRt)

# import data
loci <- rio::import("loci.tsv")
clump <- rio::import("results/clump.tsv")
meta <- rio::import("results/meta.tsv")
gwas <- rio::import("41588_2018_135_MOESM4_ESM.xlsx")
met_anno <- rio::import("41588_2018_135_MOESM7_ESM.xlsx")

# get pos37 for rsid
mart <- useMart("ENSEMBL_MART_SNP", "hsapiens_snp", host = "https://grch37.ensembl.org") # queried 
snps <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"), "snp_filter", unique(na.omit(gwas$rsID)), mart, uniqueRows = TRUE)
gwas$pos37_rsid <- snps$chrom_start[match(gwas$rsID, snps$refsnp_id)]

# prepare liftover
pos38 <- paste0("chr", gwas$Chromosome, ":", gwas$Position, "-", gwas$Position)
rio::export(as.data.frame(pos38), "stool_metabolome_hg38.txt", col.names = FALSE)

# after liftover
pos37 <- rio::import("stool_metabolome_hg19.txt", header = FALSE)
liftover_errors <- rio::import("stool_metabolome_hg19.err", format = "txt", header = FALSE)
gwas <- gwas[which(!pos38 %in% liftover_errors$V1), ]
gwas$pos37_liftover <- gsub("^.*-", "", pos37$V1)

# add pos37
gwas$Pos37 <- gwas$pos37_rsid
gwas$Pos37[which(is.na(gwas$Pos37))] <- gwas$pos37_liftover[which(is.na(gwas$Pos37))]

# clean
gwas$"Metabolite 1" <- met$Label[match(gwas$"Metabolite 1", met$Metabolite)]
gwas$"Metabolite 2" <- met$Label[match(gwas$"Metabolite 2", met$Metabolite)]
gwas$Metabolite <- gwas$"Metabolite 1"
gwas$Metabolite[which(!is.na(gwas$"Metabolite 2"))] <- paste(gwas$"Metabolite 1", gwas$"Metabolite 2", sep = " / ")[which(!is.na(gwas$"Metabolite 2"))]

split <- strsplit(meta$MarkerName, split = ":")
split <- do.call(rbind, split)
meta$chr <- split[, 1]
meta$pos <- as.numeric(split[, 2])

clump <- clump[which(clump$level == "studywide"), ]
clump$Trait <- anno$species[match(clump$trait, anno$MGS)]
clump$Trait[which(is.na(clump$Trait))] <- "Richness"

# 1 SNP per locus
order <- unique(loci$Locus)
loci <- loci[order(loci$P), ]
loci <- loci[which(!duplicated(loci$Locus)), ]
loci <- loci[match(order, loci$Locus), ]

# find stool metabolites
res <- lapply(1:nrow(loci), function(i) {
  clump <- clump[which(clump$ID == loci$rsID[i] & clump$Trait == loci$Trait[i]), ]
  snps <- c(clump$ID, unlist(strsplit(clump$SP2, split = ",")))
  chr <- loci$Chr[i]
  pos_min <- min(meta$pos[which(meta$rs_id %in% snps)])
  pos_max <- max(meta$pos[which(meta$rs_id %in% snps)])
  gwas <- gwas[which(gwas$Chromosome == chr), ]
  gwas$Locus <- loci$Locus[i]
  gwas$Trait <- loci$Trait[i]
  gwas <- gwas[order(gwas$p), ]
  gwas <- gwas[which(gwas$Position >= (pos_min - 100000) & gwas$Position <= (pos_max + 100000)), ] 
  gwas[which(!duplicated(gwas$"Metabolite")), ]
})
res <- do.call(rbind, res)
res <- res[, c("Locus", "Trait", colnames(res)[which(!colnames(res) %in% c("Locus", "Trait"))])]
res <- res[, c("Locus", "Trait", "rsID", "Chromosome", "Pos37", "Allele 0", "Allele 1", "MAF", "Metabolite", "Beta", "SE", "p")]
colnames(res) <- c("Locus", "Trait", "rsID", "Chr", "Pos37", "MajorAllele", "MinorAllele", "MAF", "Metabolite", "Beta", "SE", "P")

# export results
rio::export(res, "results/stool_metabolome.tsv")
