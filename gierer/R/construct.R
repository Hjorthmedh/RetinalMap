gierer.f <- function(rgc.dat, tec.dat,
                     nterm = 16,         #number of terminals of each RGC.
                     epsilon = 0.005, # rate of growth of r
                     kappa = 0.00000, # rate of decay of r
                     p0.fn=2         #which energy term to use, p0db
                     ) {


  ## Create the environment to return all the objects.
  g = new.env()

  g$rgc.dat = rgc.dat
  g$tec.dat = tec.dat
  ## locations of axons on retina (u,v) co-ordinates
  g$axon.positions = rgc.dat[,1:2]
  
  ## number of retinal axons
  g$axons = nrow(g$axon.positions)

  ## neighbours of axon positions
  ##axon.neighbours = neighbours.grid.8(20, 20)
  g$axon.neighbours = neighbours.voronoi(g$axon.positions, 0, 1, 0, 1)

  
  ## the axon number of each growth cone
  g$cone.axons = rep(1:g$axons, each = nterm)


  g$cones = length(g$cone.axons)  ## number of growth cones

  ## growth cones ACTIVE 1 or INACTIVE 0
  g$cone.status = rep.int(1, times = g$cones)

  ## locations of tectal positions (x,y) co-ordinates
  g$tect.positions = tec.dat[,1:2]

  ## number of tectal locations
  g$tects = nrow(g$tect.positions)

  ## neighbours of tectal positions
  g$tect.neighbours = neighbours.voronoi(g$tect.positions, 0, 1, 0, 1)

  ## tectum positions ACTIVE 1 or INACTIVE 0
  ## make sure length(which(tect.status==1)) > 1 because
  ## sample(x,1) has undesired behaviour for length(x)==1
  g$tect.status = rep.int(1, times = g$tects)

  ## which growth cone connected to which tectal position
  g$connections = sample(which(g$tect.status == 1), g$cones, replace = TRUE)

  ## initial inhibitory term
  g$r = rep.int(0, times = g$tects)

  ## The density vector
  ## This is a quick way of deriving it from the location matrix
  rho <- rep.int(0, times = g$tects)
  tL <- table(g$connections[which(g$cone.status==1)])
  names(rho) <- 1:g$tects
  rho[attributes(tL)$dimnames[[1]]] <- tL

  g$rho <- rho
  
  g$p0.fn = rep.int(p0.fn, times = g$tects)

  g$alpha = 1.0
  g$beta = 1.0
  g$sigma = 0.0
  g$tau = 0.0

  g$nterm <- nterm
  g$epsilon <- epsilon
  g$kappa <- kappa
  g$update.method = 1
  g$update.many = TRUE
  g$epoch <- 0

  g

  
}


