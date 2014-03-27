#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
stopifnot(length(args)==1)

paramfile <- args[1]
source(paramfile)

## Divert everything else now to the output file; this taken from sink() help page.
zz <- file(sprintf("%s-output.log", paramfile), open="wt")
sink(zz)
sink(zz, type="message")

## Do as much as we can after the sink, so that we can keep as many
## error messages as possible.

date()
require(WhitelawCowan)
cat(rgc.file, sep="\n")
cat(sc.file,sep="\n")
cat(op.file,sep="\n")
a <- run_WC(rgc.file,sc.file,op.file)


proc.time()
