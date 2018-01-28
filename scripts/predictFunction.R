#!/usr/bin/env Rscript
library(bmrf)
args = commandArgs(trailingOnly=TRUE)
FILEnet=args[1]
FILEann=args[2]
outprobsfile=args[3]
bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=NA, minGOsize=1,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); 
p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)