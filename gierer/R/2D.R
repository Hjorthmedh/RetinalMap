arr <- function(from, to=0) {
  if (to==0) {
    for (n in tect.neighbours[from,]) {
      if (n != 0) {
        to = n
        arrows(tect.positions[from,1], tect.positions[from,2], tect.positions[to,1], tect.positions[to,2],length=0.1)
        points(tect.positions[from,1], tect.positions[from,2], col="Blue", cex = 2, pch = 4)
        points(tect.positions[to,1], tect.positions[to,2], col="Green", cex = 2, pch = 4)
      }
    }
  } else {
    arrows(tect.positions[from,1], tect.positions[from,2], tect.positions[to,1], tect.positions[to,2])
    points(tect.positions[from,1], tect.positions[from,2], col="Blue", cex = 2, pch = 4)
    points(tect.positions[to,1], tect.positions[to,2], col="Green", cex = 2, pch = 4)
  }
}

display.locations <- function(g, main="",tectum=TRUE,
                              xmin=min(g$axon.positions[,1]), xmax=max(g$axon.positions[,1]), 
                              ymin=min(g$axon.positions[,2]), ymax=max(g$axon.positions[,2]),
                              label=FALSE,
                              ...) {
  ## show the active tectal or retinal locations

  if (tectum==TRUE) {
    tect.pos <- g$tect.positions[which(g$tect.status==1),]
    plot(tect.pos, pch=4, col="grey", cex=0.3, xlab="x", ylab="y", main=main,
         axes=F, frame.plot=T, asp=1, ...)
    ##axis(1,at=c(0,1),labels=c("Caudal","Rostral"),cex=2.5)
    axis(1,at=c(0,1),labels=c("Anterior","Posterior"),cex=2.5)
    axis(2,at=c(0,1),labels=c("Medial","Lateral"),cex=2.5)
    if (label)
      text(tect.pos, cex=0.7, labels=1:nrow(tect.pos))
  } else {
    ## Display the retina.
    axons.in.use<-unique(g$cone.axons[which(g$cone.status==1)])
    axon.pos <- g$axon.positions[axons.in.use,]
    plot(axon.pos, pch=20, cex=1, xlab="u", ylab="v", axes=F, frame.plot=T,
         xlim=c(1,0), ylim=c(1,0),      #manipulate retinal origin. TODO
         main = "The Retina",
         asp=1,
         col = colour.map(axon.pos[,1], axon.pos[,2],xmin,xmax,ymin,ymax), ...)
    axis(1,at=c(0,1),labels=c("Nasal","Temporal"),cex=2.5)
    axis(2,at=c(0,1),labels=c("Dorsal", "Ventral"),cex=2.5)
    if (label)
      text(axon.pos, cex=0.7, labels=1:nrow(axon.pos))
  }
}

## show the active axon centroids as numbers
display.numbers <- function(indices = 1:axons, points = FALSE,
                            xmin=min(axon.positions[,1]), xmax=max(axon.positions[,1]), 
                            ymin=min(axon.positions[,2]), ymax=max(axon.positions[,2])) {
  for (i in indices) {
    ## active growth cones of axon i
    in.use <- intersect(which(cone.status==1), which(cone.axons==i))
    ## where are those growth cones connected
    locs <- connections[in.use]
    ## plot
    if (points) {
      points(mean(tect.positions[locs,1]), mean(tect.positions[locs,2]), pch=16 , cex=1,
             col = colour.map(axon.positions[i,1],axon.positions[i,2],xmin,xmax,ymin,ymax))
    } else {
      text(mean(tect.positions[locs,1]), mean(tect.positions[locs,2]), labels=i, cex=1,
           col = colour.map(axon.positions[i,1],axon.positions[i,2],xmin,xmax,ymin,ymax))
    }
  }
}

## show the active "cone stars"
display.stars <- function(indices = 1:axons,
                          xmin=min(axon.positions[,1]), xmax=max(axon.positions[,1]), 
                          ymin=min(axon.positions[,2]), ymax=max(axon.positions[,2])) {
  for (i in indices) {
    ## active growth cones of axon i
    in.use <- intersect(which(cone.status==1), which(cone.axons==i))
    ## where are those growth cones connected
    locs <- connections[in.use]
    for (l in locs) {
      lines(c(tect.positions[l,1],mean(tect.positions[locs,1])), c(tect.positions[l,2], mean(tect.positions[locs,2])),
            col = colour.map(axon.positions[i,1],axon.positions[i,2],xmin,xmax,ymin,ymax), cex=0.4)
    }
  }
}

display.mesh <- function(g, indices = 1:g$axons,
                         xmin=min(g$axon.positions[,1]), xmax=max(g$axon.positions[,1]), 
                         ymin=min(g$axon.positions[,2]), ymax=max(g$axon.positions[,2])) {

  n <- g$axon.neighbours
  p <- g$axon.positions

  ignored.axons <- setdiff(1:g$axons, indices)

  n[ignored.axons,] <- 0
  n <- matrix(ifelse(n %in% ignored.axons, 0, n), nrow=nrow(n), ncol=ncol(n))

  for (i in indices) {
    ## active growth cones of axon i
    in.use <- intersect(which(g$cone.status==1), which(g$cone.axons==i))

    if (length(in.use)==0) {
      n[i,] <- 0
      n <- matrix(ifelse(n==i, 0, n), nrow=nrow(n), ncol=ncol(n))
    } else {
      ## where are those growth cones connected
      locs <- g$connections[in.use]
      
      p[i,1] <- mean(g$tect.positions[locs,1])
      p[i,2] <- mean(g$tect.positions[locs,2])
    }
  }

  for (r in 1:nrow(n)) {
    for (c in which(n[r,]!=0)) {
      lines( c(p[r,1] , p[n[r,c],1]) , c(p[r,2] , p[n[r,c],2]) ,
            col = colour.map(g$axon.positions[r,1],g$axon.positions[r,2],
              xmin,xmax,ymin,ymax) )
    }
  }
  ## After drawing the mesh, keep track of the centre-of-mass for each RGC target.
  g$com <- p
}

display.potential <- function(axon) {
  require(lattice)

  u <- axon.pot.positions[axon,1]
  v <- axon.pot.positions[axon,2]
  x <- tect.pot.positions[,1]
  y <- tect.pot.positions[,2]

  g <- cbind(x=tect.positions[,1], y=tect.positions[,2])
  g <- as.data.frame(g)
  g$z <- p.gierer2(x,y,u,v)
  g$gr <- 1

  wireframe(z ~ x * y, data = g, groups = gr, drape=T, scales = list(arrows = FALSE), cex=0.4, zlab="p", col.regions = terrain.colors(100), main="alpha = 100")
}

p.gierer2 <- function(x, y, u, v) {
  (exp(-alpha*u) / exp(-alpha*x)) + (exp(-alpha*x) / exp(-alpha*u)) +
    (exp(-beta*v) / exp(-beta*y)) + (exp(-beta*y) / exp(-beta*v)) 
}

display.r <- function() {
  require(rgl)

  plot3d(tect.positions[,1], tect.positions[,2], r, xlab = "x", ylab = "y", zlab = "r")
}

display.rho <- function() {
  require(rgl)

  plot3d(tect.positions[,1], tect.positions[,2], rho, xlab = "x", ylab = "y", zlab = "rho")
}

## colour map projects 2D x-y plane to plane cross section of RGB cube
## one of the 3 planes that avoids the (0,0,0) and (255,255,255) corners
colour.map <- function(x, y, xmin, xmax, ymin, ymax) {
  ## map x and y values to [0,1] interval
  x <- (x - xmin) / (xmax - xmin)
  if (ymin == ymax) {	
    y <- 0
  } else {
    y <- (y - ymin) / (ymax - ymin)
  }
  rgb(1-x, x, y)
}

iterate <- function(g, iterations, sequence = rep.int(0,times=sum(g$cone.status))) {
  output <- 
    .C("iteratedb",
               as.integer(g$axons),
               as.double(g$rgc.dat),
               as.integer(g$cones),
               as.integer(g$cone.axons),
               as.integer(which(g$cone.status==1)),
               as.integer(sum(g$cone.status)),
               as.integer(g$tects),
               as.double(g$tec.dat),
               as.integer(ncol(g$tect.neighbours)),
               as.integer(g$tect.neighbours),
               as.integer(g$tect.status),
               connections.out = as.integer(g$connections),
               r.out = as.double(g$r),
               as.double(g$epsilon),
               as.double(g$kappa),
               rho.out = as.integer(g$rho),
               as.integer(g$p0.fn),
               as.double(g$alpha),
               as.double(g$beta),
               as.double(g$sigma),
               as.double(g$tau),
               as.integer(iterations),
               as.integer(sequence),
               as.integer(length(sequence)),
               as.integer(g$update.method),
               as.integer(g$update.many))

  g$connections <- output[["connections.out"]]
  g$r <- output[["r.out"]]
  g$rho <- output[["rho.out"]]
  g$epoch <- g$epoch + iterations
}

## Inactivate the growth cones of particular axons
kill.axons <- function(indices) {
  for (i in indices) {
    cones.of.axon <- which(cone.axons == i)
    kill.cones(cones.of.axon)
  }
}

## Activate the growth cones of particular axons
## By default, the revived cones are spread randomly across the whole tectum
revive.axons <- function(indices) {
  for (i in indices) {
    cones.of.axon <- which(cone.axons == i)
    revive.cones(cones.of.axon, "tectum")
  }
}

## Inactivate particular growth cones
kill.cones <- function(indices) {
  kill.list <- intersect(which(cone.status == 1), indices)
  cone.status[kill.list] <<- 0

  ## The density vector
  ## This is a quick way of deriving it from the location matrix
  rho <<- rep.int(0, times = tects)	
  tL <- table(connections[which(cone.status==1)])
  names(rho) <<- 1:tects
  rho[attributes(tL)$dimnames[[1]]] <<- tL
}

## Activate particular growth cones
## Option to distribute growth cones randomly over active tectum or randomly across pre-connected cone locations
revive.cones <- function(indices, option = "tectum") {
  revive.list <- intersect(which(cone.status == 0), indices)	

  ## Distribute over entire active tectum
  if (option == "tectum") {
    cone.status[revive.list] <<- 1		
    connections[revive.list] <<- sample(which(tect.status == 1), length(revive.list), replace = TRUE)
  }
  ## Distribute over pre-connected locations
  if (option == "axon") {
    for (c in revive.list) {
      axon.number <- cone.axons[c]
      active.gcs <- intersect(which(cone.status == 1), which(cone.axons == axon.number))
      active.gcs.connections <- connections[active.gcs]

      ## No other connections, then do over entire tectum
      if (length(active.gcs.connections) == 0) {
        connections[c] <<- sample(which(tect.status == 1), 1)
        ## Otherwise randomly over current connections
      } else {
        ## sample(x,1) has undesired behaviour for length(x)==1
        if (length(active.gcs.connections) == 1) {
          connections[c] <<- active.gcs.connections
        } else {
          connections[c] <<- sample(active.gcs.connections, 1)
        }
      }
    }
    cone.status[revive.list] <<- 1		
  }
}

## Inactivate particular tectal locations
## kill.whole.axon = TRUE kills all cones of the axons affected
kill.tects <- function(indices, kill.whole.axon = FALSE) {
  kill.list <- intersect(which(tect.status == 1), indices)
  tect.status[kill.list] <<- 0
  
  if (kill.whole.axon == FALSE) {
    kill.cones(which((connections %in% kill.list) == TRUE))
  }
  if (kill.whole.axon == TRUE) {
    kill.axons(unique(cone.axons[which((connections %in% kill.list) == TRUE)]))
  }
}

## Activate particular tectal locations
revive.tects <- function(indices) {
  revive.list <- intersect(which(tect.status == 0), indices)
  tect.status[revive.list] <<- 1
}

## returns the clockwise rotation of points found in the given window
rotate <- function(pot.positions, xmin, xmax, ymin, ymax, rads) {
  centre.point <- c(mean(c(xmax,xmin)), mean(c(ymax,ymin)))

  ## get the points that are within specified window
  to.rotate <- intersect(		intersect(which(pot.positions[,1] >= xmin), which(pot.positions[,1] <= xmax)) ,
                         intersect(which(pot.positions[,2] >= ymin), which(pot.positions[,2] <= ymax)) )

  ## for clockwise rotation
  rotation.matrix <- matrix(c(cos(rads),sin(rads),-sin(rads),cos(rads)), nrow=2, ncol=2, byrow=TRUE)
  
  ## map points to centre.point as origin
  translated.points <- t(pot.positions[to.rotate,]) - centre.point

  ## do rotation
  rotated.points <- rotation.matrix %*% translated.points

  ## map points back to normal locations
  translated.points <- t(rotated.points + centre.point)

  ## update the affected values
  pot.positions[to.rotate,] <- translated.points

  pot.positions
}


######################################################################
## These functions moved from Environment.R

### FUNCTIONS ###

## use to construct rectangular grid
## leave out ymax for regular grid
construct.grid <- function(rows, cols, xmin, xmax, ymin,
                           ymax = (((xmax-xmin)/(cols-1))*(rows-1) + ymin)) {
  ## initialise
  co.ords <- matrix(0, nrow = rows*cols, ncol = 2)

  ## x-ordinates
  co.ords[,1] = rep(seq(from = xmin, to = xmax, length = cols), times = rows)

  ## y-ordinates
  co.ords[,2] = rep(seq(from = ymin, to = ymax, length = rows), each = cols)

  ## return
  co.ords
}

## use to construct hexagonal grid
## leave out ymax for regular grid
construct.hex <- function(rows, cols, xmin, xmax, ymin,
                          ymax = (((xmax-xmin)/(cols-1))*(rows-1) + ymin)*(sqrt(3)/2)) {
  ## initialise
  co.ords <- matrix(0, nrow = rows*cols, ncol = 2)

  ## find gap length in x direction
  x.gap <- (xmax-xmin)/(cols-1)

  ## x-ordinates
  to.add <- ((((1:(rows*cols)) - 1) %/% cols) %% 2)*(x.gap/2)
  co.ords[,1] = rep(seq(from = xmin, to = xmax-x.gap/2, length = cols), times = rows) + to.add

  ## y-ordinates
  co.ords[,2] = rep(seq(from = ymin, to = ymax, length = rows), each = cols)

  ## return
  co.ords
}

## use to construct random uniform points
construct.random <- function(N, xmin, xmax, ymin, ymax) {
  ## initialise
  co.ords <- matrix(0, nrow = N, ncol = 2)

  ## x-ordinates
  co.ords[,1] = runif(N, xmin, xmax)

  ## y-ordinates
  co.ords[,2] = runif(N, ymin, ymax)

  ## return
  co.ords
}

## use to construct dmin points
construct.dmin <- function(N, xmin, xmax, ymin, ymax, mean, sd) {
  require(sjedmin)
  window <- c(xmin, xmax, ymin, ymax)
  d = dminlul(window, N, mean, sd)
  cbind(d$x, d$y)
}

## use for rectangular grid and 4 neighbours
neighbours.grid.4 <- function(rows, cols) {
  total.tects <- rows*cols

  ## initialise
  neighbour.indices <- matrix(0, nrow = total.tects, ncol = 4)

  for (i in 1:total.tects) {
    ## not on bottom row
    if (i > cols) {neighbour.indices[i,1] = i - cols}
    ## not on top row
    if (i <= total.tects - cols) {neighbour.indices[i,2] = i + cols}
    ## not on left col
    if ((i %% cols) != 1) {neighbour.indices[i,3] = i - 1}
    ## not on right col
    if ((i %% cols) != 0) {neighbour.indices[i,4] = i + 1}
    
    neighbour.indices[i,] = sort(neighbour.indices[i,], decreasing = TRUE)
  }

  ## return
  neighbour.indices
}

## use for rectangular grid and 8 neighbours
neighbours.grid.8 <- function(rows, cols) {
  total.tects <- rows*cols

  ## initialise
  neighbour.indices <- matrix(0, nrow = total.tects, ncol = 8)

  for (i in 1:total.tects) {
    ## not on bottom row
    if (i > cols) {neighbour.indices[i,1] = i - cols}
    ## not on top row
    if (i <= total.tects - cols) {neighbour.indices[i,2] = i + cols}
    ## not on left col
    if ((i %% cols) != 1) {neighbour.indices[i,3] = i - 1}
    ## not on right col
    if ((i %% cols) != 0) {neighbour.indices[i,4] = i + 1}
    ## not on bottom row or left col	
    if ((i > cols)&&((i %% cols) != 1)) {neighbour.indices[i,5] = i - cols - 1}
    ## not on bottom row or right col
    if ((i > cols)&&((i %% cols) != 0)) {neighbour.indices[i,6] = i - cols + 1}
    ## not on top row or left col
    if ((i <= total.tects - cols)&&((i %% cols) != 1)) {neighbour.indices[i,7] = i + cols - 1}
    ## not on top row or right col
    if ((i <= total.tects - cols)&&((i %% cols) != 0)) {neighbour.indices[i,8] = i + cols + 1}

    neighbour.indices[i,] = sort(neighbour.indices[i,], decreasing = TRUE)
  }

  ## return
  neighbour.indices
}

## use for hexagonal grid and 6 neighbours
neighbours.hex.6 <- function(rows, cols) {
  total.tects <- rows*cols

  ## initialise
  neighbour.indices <- matrix(0, nrow = total.tects, ncol = 6)

  for (i in 1:total.tects) {
    ## not on bottom row
    if (i > cols) {neighbour.indices[i,1] = i - cols}
    ## not on top row
    if (i <= total.tects - cols) {neighbour.indices[i,2] = i + cols}
    ## not on left col
    if ((i %% cols) != 1) {neighbour.indices[i,3] = i - 1}
    ## not on right col
    if ((i %% cols) != 0) {neighbour.indices[i,4] = i + 1}

    offset.row <- ((i-1) %/% cols) %% 2

    if (offset.row == 0) {
      ## not on bottom row or left col	
      if ((i > cols)&&((i %% cols) != 1)) {neighbour.indices[i,5] = i - cols - 1}
      ## not on top row or left col
      if ((i < total.tects - cols + 1)&&((i %% cols) != 1)) {neighbour.indices[i,6] = i + cols - 1}
    } else {
      ## not on bottom row or right col
      if ((i > cols)&&((i %% cols) != 0)) {neighbour.indices[i,5] = i - cols + 1}
      ## not on top row or right col
      if ((i < total.tects - cols + 1)&&((i %% cols) != 0)) {neighbour.indices[i,6] = i + cols + 1}
    }

    neighbour.indices[i,] = sort(neighbour.indices[i,], decreasing = TRUE)
  }

  ## return
  neighbour.indices
}

## use for calculation of neighbours based on euclidean distance
## positions as returned by the construct functions above
## num.nearest is for integer number of neighbours specified
## threshold is for maximum distance that a neighbour can be
neighbours.nearest <- function(positions, num.nearest = 0, threshold = 0) {
  total <- nrow(positions)

  ## select all as neighbours
  if ((num.nearest == 0)&&(threshold == 0)) {
    ## initialise
    neighbour.indices <- matrix(0, nrow = total, ncol = total - 1)

    for (r in 1:total) {
      ## all are neighbours except r			
      neighbour.indices[r,] <- (1:total)[-r]
    }

  }
  
  ## all within a threshold
  if ((num.nearest == 0)&&(threshold != 0)) {
    ## find max size of neighbours matrix
    max.size <- 0
    for (r in 1:total) {
      distances <- sqrt(colSums((t(positions) - positions[r,])^2))
      num.found <- length(which(distances <= threshold)) - 1
      if (num.found > max.size) {max.size <- num.found}	
    }

    ## initialise
    neighbour.indices <- matrix(0, nrow = total, ncol = max.size)	
    
    if (max.size == 0) {break}
    
    for (r in 1:total) {
      distances <- sqrt(colSums((t(positions[-r,]) - positions[r,])^2))
      found <- which(distances <= threshold)
      neighbour.indices[r,1:(length(found))] <- ifelse(found >= r, found + 1, found)
    }
  }

  ## no threshold
  if ((num.nearest != 0)&&(threshold == 0)) {
    ## initialise
    neighbour.indices <- matrix(0, nrow = total, ncol = num.nearest)	

    for (r in 1:total) {
      distances <- sqrt(colSums((t(positions) - positions[r,])^2))
      neighbour.indices[r,] <- sort(distances, index.return = TRUE)$ix[2:(num.nearest+1)]
    }
  }

  ## threshold first, then num nearest
  if ((num.nearest != 0)&&(threshold != 0)) {
    ## initialise
    neighbour.indices <- matrix(0, nrow = total, ncol = num.nearest)	

    for (r in 1:total) {
      distances <- sqrt(colSums((t(positions) - positions[r,])^2))
      found <- (which(distances <= threshold))
      neighbour.indices[r,] <- found[sort(distances[found], index.return = TRUE)$ix[2:(num.nearest+1)]]
    }
  }

  ## replace NAs with 0
  neighbour.indices[is.na(neighbour.indices)] <- 0

  ## return
  neighbour.indices
}

neighbours.voronoi <- function(positions, xmin, xmax, ymin, ymax) {
  require(sjevor)
  window <- c(xmin, xmax, ymin, ymax)
  ## include 'rejects' by using opts="asi"
  v <- vorcr2(positions, window, opts="asi")
  ifelse(v$neighs == -1, 0, v$neighs)
}


write.connections <- function(g, opfile) {
  conn <- with(g, 
    sapply(1:axons, function(i) tabulate(connections[which(cone.axons ==i)], nbins=tects)))
  write.table(conn, file=opfile, sep=',', row.names=F, col.names=F)
  stopifnot(isTRUE(all.equal(colSums(conn), rep(16, g$axons))))
}

write.targets <- function(g, opfile) {
  ## Similar to write.connections, but just output one row per RGC listing
  ## the targets of that RGC.  If the RGC makes contact with the same tectal
  ## cell multiple times, it is listed multiple times.
  conn <- t(with(g, sapply(1:axons, function(i) 
                         connections[which(cone.axons == i)])))
  
  stopifnot(ncol(conn) == g$nterm)
  write.table(conn, file=opfile, sep=',', row.names=F, col.names=F)
}

display.map <- function(g, main=NULL, label=FALSE) {

  if(is.null(main)) {
    main <- sprintf("epoch %d", g$epoch)
  }

  display.locations(g, main=main, label=label)
  display.mesh(g)
}

