#!/usr/bin/env Rscript
# This script expects one filename to be passed, the parameter file.
# This parameter file (which is just a R script) must stores the values of 
# 3 variables:
# 
# rgc.file
# sc.file
# op.file
#
# Other parameters (e.g. kappa, epsilon) can be set.
#

## The name of the paramfile is used as the prefix for output files,
## including a log file,  and a pdf showing how the map develops.
args <- commandArgs(TRUE)
stopifnot(length(args)==1)

## Defaults for key parameters; these can be over-ridden in "paramfile".
p0.fn <- 4
kappa <- 0.0
epsilon <- 0.005

## Read in the values from the parameter file.
paramfile <- args[1]
source(paramfile)

## Divert everything else now to the output file; this taken from
## sink() help page. The log file is useful for seeing how long the
## simulation is taking.
zz <- file(sprintf("%s-output.log", paramfile), open="wt")
sink(zz)
sink(zz, type="message")

## Do as much as we can after the sink, so that we can keep as many
## error messages as possible.

date()
require(gierer)
rgc.dat <- as.matrix(read.table(rgc.file))
tec.dat <- as.matrix(read.table(sc.file))

a <- gierer.f(rgc.dat, tec.dat, epsilon=epsilon, kappa=kappa, p0.fn=p0.fn)

summary.gierer(a)


pdf(file=sprintf("%s.pdf", paramfile), width=11, height=8)
on.exit(dev.off())

par(mfrow=c(2,3))
display.map(a); showxyz(a, "rho"); showxyz(a, "r")
com <- a$com
epoch <- 0
##checkpoints <- c(0, 5, 20, 50, 100, 500, 1000, 1100, 1200, 1300, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000)
##checkpoints <- seq(from=0, by=100, to=1500)
#checkpoints <- c(0, 100)
checkpoints <- c(0, 5, 20, 100, 500, 1000, 3000, 5000, 10000)
num.checkpoints <- length(checkpoints)
com.dist <- rep(0, num.checkpoints)
for (i in 2:num.checkpoints) {
  cat(sprintf("Epoch %d %s\n", epoch, date()))
  epochs <- checkpoints[i] - checkpoints[i-1]
  iterate(a, epochs)
  epoch <- epoch + epochs
  main.label <- sprintf("Ep %d pfn %d\n", epoch, a$p0.fn[1])
  display.map(a,main=main.label); showxyz(a, 'rho'); showxyz(a, 'r')
  com.dist[i] <- com.distance(com, a$com)
  com <- a$com
}

plot(checkpoints, com.dist,
     ylab='Mean RGC CoM movement since last checkpoint', type='b')
display.locations(a,tectum=FALSE, label=FALSE)

## do you want a full (connections) or sparse (targets) matrix?
##write.connections(a, op.file)
write.targets(a, op.file)

date()
proc.time()
