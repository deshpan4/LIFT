library(bmrf)
FILEnet="L-P-P-P-matrix.txt"
FILEann="P-GOTERM-BMRFFILE.txt"
outprobsfile="lncRNAs-function-pred-DGE.txt"
bmrf <- read.bmrf(net.file=FILEnet, go.file=FILEann, fd.file=NA, minGOsize=1,minFDsize=3, maxGOsize=Inf, maxFDsize=Inf )
set.seed(10); p <- predict.bmrf(bmrf, file=outprobsfile, format="3col", verbose=TRUE)
