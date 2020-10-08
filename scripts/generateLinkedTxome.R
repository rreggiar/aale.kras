#!/user/bin/env Rscript

# rreggiar@ucsc.edu

# [guide to generating with salmon](https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/)

suppressPackageStartupMessages({
  library(tximeta)
})

paths.in <- paths.in <- scan(file=file("stdin", "r"), what="character", n=1)

if (length(paths.in) < 1) {
  cat('please provide all arguments: [1] input directory')
  quit(status = 1)
}

input.dir <- paths.in[1]
gencode.version <- substring(input.dir, first = nchar(input.dir)-1, last = nchar(input.dir))

tximeta::makeLinkedTxome(
    indexDir = file.path(input.dir, paste0('sel.align.gencode.v', gencode.version, '.process.aware.salmon.v', '*.sidx')),
    source = "GENCODE", genome = "GRCh38",
    organism = "Homo sapiens", release = "Hg38",
    fasta = file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.expanded.fa')),
    gtf = file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.expanded.gtf')),
    write = TRUE, jsonFile = file.path(input.dir, paste0('gencode.v', gencode.version, '.annotation.expanded.json'))
)