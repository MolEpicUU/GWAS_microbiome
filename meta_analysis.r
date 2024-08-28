# load libraries
library(rio)

# prepare METAL input files
files <- list.files("data/scapis/")
files <- gsub(".tsv.gz", "", files)
lapply(files, function(x) {
  text <- c(
    "SEPARATOR COMMA",
    "MARKER ID",
    "ALLELE ALLELE1 ALLELE0",
    "FREQ A1FREQ",
    "EFFECT BETA",
    "STDERR SE",
    "SCHEME STDERR",
    "AVERAGEFREQ ON",
    "MINMAXFREQ ON",
    paste0("PROCESS data/scapis/", x, ".tsv.gz"),
    paste0("PROCESS data/cosmc/", x, ".tsv.gz"),
    paste0("PROCESS data/smcc/", x, ".tsv.gz"),
    paste0("PROCESS data/mos/", x, ".tsv.gz"),
    paste0("OUTFILE output/", x, "_ .txt"),
    "ANALYZE HETEROGENEITY"
  )
  writeLines(text, paste0("input/", x, ".txt"))
})

# run METAL
files <- list.files("input", full.names = TRUE)
chunks <- split(files, cut(seq_along(files), 20, labels = FALSE))
lapply(seq_along(chunks), function(i) {
  lapply(chunks[[i]], function(x) {
    system(paste0("sbatch -t 02:00:00 -J mtl_", i, " -d singleton -o logs/metal_", gsub("input/", "", x), '.log --wrap="ml bioinfo-tools; ml METAL; metal ', x, '"'))
  })
})
