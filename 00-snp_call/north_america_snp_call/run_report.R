#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
if (length(args) < 2) stop("Bad args, usage refdir cmpdir")

rmd_file  <- args[1]
gt        <- args[2]
ao        <- args[3]
ro        <- args[4]
qual      <- args[5]
scaff_len <- args[6]
q         <- args[7]
samples   <- args[8]

out_dir   <- paste(getwd(), args[9], sep="/")
fig_dir   <- "figs"

write(paste("Wrapper is ran on:", getwd()),
      stderr())

library(ezknitr)
ezknit(file    = rmd_file,
       out_dir = out_dir,
       fig_dir = fig_dir,
       params  = list("gt"          = gt,
                      "ao"          = ao,
                      "ro"          = ro,
                      "qual"        = qual,
                      "q"           = q,
                      "samples"     = samples,
                      "scaff_len"   = scaff_len),
       verbose = TRUE)
