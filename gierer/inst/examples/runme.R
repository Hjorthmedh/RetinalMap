require(gierer)
rgc.file <- system.file('examples/Gierer2D-500-WT-retinal-gradients.txt', package='gierer')
sc.file <- system.file('examples/Gierer2D-500-WT-SC-gradients.txt', package='gierer')
op.file <- "/tmp/wts.dat"

rgc.dat <- as.matrix(read.table(rgc.file))
tec.dat <- as.matrix(read.table(sc.file))

a <- gierer.f(rgc.dat, tec.dat)

par(mfrow=c(2,3))
display.locations(a,tectum=FALSE, label=FALSE)
display.map(a)

epoch <- 0
##checkpoints <- c(0, 5, 20, 50, 100, 500, 1000, 2000)
checkpoints <- c(0, 100)
checkpoints <- c(0, 5, 20, 50, 100)
num.checkpoints <- length(checkpoints)
for (i in 2:num.checkpoints) {
  cat(sprintf("Epoch %d\n", epoch))
  epochs <- checkpoints[i] - checkpoints[i-1]
  iterate(a, epochs)
  epoch <- epoch + epochs
  display.map(a, label=FALSE)
}

## do you want a full (connections) or sparse (targets) matrix?
##write.connections(a, op.file)
write.targets(a, op.file)
