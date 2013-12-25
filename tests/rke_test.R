##   Testing for rke.* functions.

test.rke.extended.Rsubgraph <- function() {
  #
  # Tests for the extended R subgraph function.
  #
  pool = rrke.pool(m=3, n=30, uniform.pra=T)
  legit.pairs = c("O-A", "A-O", "O-B", "B-O", "A-B", "B-A")
  rke = pool$rke.all
  xRKE = rke.extended.Rsubgraph(rke)
  CHECK_MEMBER(xRKE$pairs$desc, legit.pairs)
  print(sprintf("Total edges %d", nrow(xRKE$edges)))
  # Check all OU pairs -> UD pairs of the same hospital
  for(ou in c("O-A", "O-B")) {
    pairs = subset(xRKE$pairs, desc==ou)
    for(id in pairs$pair.id) {
      h1 = subset(xRKE$pairs, pair.id==id)$hospital
      to.pair.ids = subset(xRKE$edges, pair.id1==id & can.donate==1)$pair.id2
      to.pairs.desc = subset(xRKE$pairs, pair.id %in% to.pair.ids)$desc
      # 1: O-A -> B-O  only so that a virtual pair is created.
      if(ou=="O-A") CHECK_TRUE(all(to.pairs.desc=="B-O"))
      else CHECK_TRUE(all(to.pairs.desc=="A-O"))
      h2 = subset(xRKE$pairs, pair.id %in% to.pair.ids)$hospital
      
      # 2: OD -> UD should only be within hospital
      CHECK_SETEQ(h2, h1)
      # 3: Outlinks in the extended graph is subset of the outlinks in the original graph
      #    (i.e make sure we don't introduce new links)
      to.pair.ids.original = subset(rke$edges, pair.id1==id & can.donate==1)$pair.id2
      CHECK_MEMBER(to.pair.ids, to.pair.ids.original)
    }
  }
  # CHECK that UD pairs  -> (A-B or B-A) 
  for(ud in c("A-O", "B-O")) {
    pairs = subset(xRKE$pairs, desc==ud)
    for(id in pairs$pair.id) {
      to.pair.ids = subset(xRKE$edges, pair.id1==id & can.donate==1)$pair.id2
      to.pairs.desc = subset(xRKE$pairs, pair.id %in% to.pair.ids)$desc
      # 2 : A-O -> B-A  only (any hospital).
      if(ud=="A-O") CHECK_TRUE(all(to.pairs.desc=="B-A"))
      else CHECK_TRUE(all(to.pairs.desc=="A-B"))
      
      # 3: Outlinks in the extended graph is subset of the outlinks in the original graph
      #    (i.e make sure we don't introduce new links)
      to.pair.ids.original = subset(rke$edges, pair.id1==id & can.donate==1)$pair.id2
      CHECK_MEMBER(to.pair.ids, to.pair.ids.original)
    }
  }
}

get.RKE0 <- function(ntriangles) {
  # Creates an RKE with a random #triangles
  # If N = #triangles then we know max matching = 3 * N
  # If only 2-way exchanges are used, then we know this is 2 * N
  # See a picture at "tests/RKE0.png"
  if (ntriangles == 0)
    return(empty.rke())
  pairs = kPairs[sample(1:nrow(kPairs), size=3 * ntriangles, replace=T),]
  pairs$pra <- rep(0, nrow(pairs))
  pairs$pair.id <- 1:nrow(pairs)
  pairs$hospital <- rep(1, nrow(pairs))
  rke = empty.rke()
  rke$pairs = pairs
  rke <- rke.update.new.pairs(rke, keep.edges=c())
  filter <- function(x) {
    i = intersect(which(x < 3 * ntriangles), which(x >= 1))
    return(x[i])
  }
  give.to <- function(i) {
    rem = i %% 3
    if(rem == 1) return(filter(c(i+1, i-3, i+3)))
    if(rem == 2) return(c(i+1))
    if(rem == 0) return(c(i-2))
  }
  rke$edges$can.donate <- rep(0, nrow(rke$edges))
  for(i in rke.pair.ids(rke)) {
    j = intersect(which(rke$edges$pair.id1 == i), which(rke$edges$pair.id2 %in% give.to(i)))
    rke$edges$can.donate[j] <- 1
  }
  rke$A = map.edges.adjacency(edges=rke$edges, all.pair.ids=rke.pair.ids(rke))
  return(rke)
}

test.rke.operations <- function() {
  npairs = rpois(1, lambda=80)
  rke = rrke(npairs)
  loginfo("Checking RKE")
  CHECK_rke(rke)
  # rke.size, rke.remove.pairs, rke.keep.pairs, rke.add, rke.subgraph,
  # rke.edge.ids, rke.pair.ids, rke.hospital.ids, rke.pairs.hospitals,
  # rke.filter.pairs, pairid.to.pc, pairid.to.hospital
  loginfo("Checking rke.size()")
  CHECK_EQ(rke.size(rke), npairs)  # rke.size
  loginfo("Checking rke.pair.ids()")
  CHECK_EQ(rke.pair.ids(rke=rke), 1:npairs)  # rke.pair.ids
  rm.pairs = sample(1:npairs, as.integer(0.2 * npairs), replace=F)
  rke1 <- rke.remove.pairs(rke, rm.pairs)  # rke.remove.pairs
  CHECK_EQ(rke.size(rke1), npairs-length(rm.pairs))
  CHECK_SETEQ(union(rke.pair.ids(rke1), rm.pairs), 1:npairs)
  loginfo("Check whether removed RKE is subset")
  CHECK_rke_subset(rke.smaller=rke1, rke.bigger=rke)
  loginfo("Checking rke.subgraph()")
  # pair type should be one of {O, U, R, S}
  CHECK_ERROR(rke.subgraph(rke, pair.type="NONE"))
  x = rke.subgraph(rke, pair.type="U")
  CHECK_EQ(rke.size(x), nrow(subset(rke$pairs, pair.type=="U")))
  CHECK_rke_subset(rke.smaller=x, rke.bigger=rke)
  loginfo("Testing pairid.to.pc()")
  sample.ids <- sample(rke$pairs$pair.id, size=rpois(1, lambda=2), replace=F)
  x <- subset(rke$pairs, pair.id %in% sample.ids)
  CHECK_SETEQ(x$pc, pairid.to.pc(rke, sample.ids))
  CHECK_SETEQ(x$hospital, pairid.to.hospital(rke, sample.ids))
  
  loginfo("Checking filter() functions")
  CHECK_ERROR(rke.filter.pairs(rke, attrs=c("asas"), values=c(10)))
  CHECK_ERROR(rke.filter.pairs(rke, attrs=c("pc"), values=c()))
  random.pc <- sample(kPairCodes, 1)
  x = rke.filter.pairs(rke, attrs=c("pc"), values=c(random.pc))
  y = subset(rke$pairs, pair.id %in% x)$pc
  if(length(y) > 0) CHECK_SETEQ(y, c(random.pc))
  
  loginfo("Checking rke.add()...")
  
  npairs <- 250
  rke1 = rrke(npairs, uniform.pra=F)
  rke2 = rrke(npairs, uniform.pra=F, pair.ids=(npairs+1:npairs))
  ## Pick a pair at random.
  random.pair <- sample(1:nrow(rke1$pairs), 1)
  p1 <- rke1$pairs[random.pair,]$pair.id
  incoming.pcs <- subset(kPairs, rke1$pairs[random.pair,]$patient %% donor == 0)$pc
  # find all pair ids in RKE 2 that can give to the random pair
  p2 <- subset(rke2$pairs, pc %in% incoming.pcs)$pair.id
  pra2 <- subset(rke2$pairs, pair.id %in% p2)$pra
  random.pra <- rke1$pairs[random.pair,]$pra
  loginfo(sprintf("Sampled PRA %.4f", random.pra))
  edge.probs <- rep(1-random.pra, length(p2))
  CHECK_EQ(length(edge.probs), length(pra2))
  loginfo(sprintf("Found %d matching pairs ", length(p2)))
  rke = rke.add(rke1, rke2)
  final.edges <- subset(rke$edges, pair.id2 == p1 & pair.id1 %in% p2)$can.donate
  xmean <- mean(final.edges)
  vtheor <- sum(edge.probs * (1-edge.probs)) / length(edge.probs)^2
  loginfo(sprintf("Observed mean=%.4f Theor=%.4f  V=%.3f", xmean, mean(edge.probs), vtheor))
  CHECK_NEAR(xmean, mean(edge.probs), tol=3 * sqrt(vtheor))
  
  loginfo("Transform everyone to be very sensititive")
  rke1$pairs$pra <- 1
  rke2$pairs$pra <- 1
  # since PRA=0 for everyone, there will be no extra edges
  get.no.edges <- function(rke) {
    sum(rke$edges$can.donate)
  }
  rke.both <- rke.add(rke1, rke2)
  CHECK_rke_subset(rke1, rke.both)
  CHECK_rke_subset(rke2, rke.both)
  CHECK_EQ(get.no.edges(rke.both), get.no.edges(rke1) + get.no.edges(rke2))
}

test.rke.cycles <- function() {
  ntriangles = rpois(1, lambda=20)
  rke = get.RKE0(ntriangles)
  loginfo(sprintf("Testing with %d triangles..", ntriangles))
  cycles2 = rke.cycles(rke, include.3way=F)
  cycles3 = rke.cycles(rke, include.3way=T)
  CHECK_EQ(nrow(cycles2), ntriangles-1)
  CHECK_EQ(nrow(subset(cycles3, type==3)), ntriangles)
}
