#!/usr/bin/env Rscript
library(scarHRD)
args = commandArgs(trailingOnly=TRUE)

# run scarHRD
scarHRD_input <- args
scar_score(scarHRD_input, reference = "grch38", seqz=FALSE, chr.in.names=TRUE)

