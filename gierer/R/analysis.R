com.distance <- function(old, cur) {
  ## Compare old and new COM distance.
  delta <- old - cur
  dists <- apply(delta, 1, function(x) {sqrt(sum(x**2))})
  mean(dists)
}

scatter.grad <- function(x, y, z, lo, hi, ncols=20, main, leg.title='title', ...) {
  ## Plot 'z' at irregular x,y locations by circles of varying grey levels.
  if(missing(lo)) lo <- min(z)
  if(missing(hi)) hi <- max(z)

  ## check when all values are exactly the same, e.g. rep(0, 100)
  if( hi-lo < 1e-6)
    hi <- lo+0.001
  z.cut <- cut(z, breaks=seq(from=lo, to=hi, length=ncols), include.lowest=TRUE)
  cols=gray(seq(from=0, to=1, length=ncols))
  plot(x, y, pch=21, bg=cols[z.cut], asp=1, lwd=0.2, ...)
  legend('topleft',legend=signif(c(hi, (hi+lo)/2, lo),5), pch=21,
         pt.bg=cols[c(ncols, ncols/2,1)])
  if (missing(main)) {
    title <- sprintf('%s: mean %.3f sd %.3f\n', leg.title, mean(z), sd(z))
    title(main=title)
  } else {
    title(main=main)
  }
}

showxyz <- function(a, name, title=name) {
  ## if(missing(title))
  ##   title <- name
  scatter.grad(a$tect.positions[,1], a$tect.positions[,2], a[[name]],
               leg.title=title)}


summary.gierer <- function(a) {
  ## Short summary of the object.
  cat(sprintf("%d -> %d neurons\n", nrow(a$rgc.dat), nrow(a$tec.dat)))
  cat(sprintf("p0.fn = %d\n", a$p0.fn[1]))
  cat(sprintf("epsilion = %f\n", a$epsilon))
  cat(sprintf("kappa = %f\n", a$kappa))
}
