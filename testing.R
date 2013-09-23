# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
source("../r-toolkit/checks.R")

repeat.test = function(test, ntrials) {
  for(t in 1:ntrials) {
    res = test()
    print(sprintf("t=%3d/%3d:   [OK]", t, ntrials))
  }  
}

# Testing of terminology.R
test.terminology.blood <- function() {
  # Tests functions of blood code and types.
  CHECK_EQ(as.blood.code(c("O", "O", "A", "AB", "B", "B")), c(1,1,2,6,3,3))
  CHECK_ERROR(get.blood.code.prob(blood.code=5), msg="No blood code 5 exists.")
  CHECK_EQ(as.blood.code(rep(kBloodTypes, 3)), rep(kBloodCodes, 3))
  print("[PASS]")
}

####    TESTS for terminology.R
test.terminology.kpairs <- function() {
  # kPairs = all 16 pair code combinations
  # Here we test whether the information in this data frame is correct.
  #
  # TEST #pairs = #code x #code
  CHECK_EQ(nrow(kPairs), length(kBloodCodes)^2)
  # Pick a pair code at random
  pcode = sample(kPairCodes, 1)
  pc.combinations = expand.grid(patient=kBloodCodes, donor=kBloodCodes)
  # TEST whether (donor, patient) | pc   are correct.
  CHECK_EQ(subset(kPairs, pc==pcode, select=c(donor, patient)),
           pc.combinations[pcode, ])
  x = subset(kPairs, pc==pcode)
  # Check whether the information in kPairs is correct for this particular
  # pair code.
  donor.id = which(kBloodCodes==x$donor)
  patient.id = which(kBloodCodes==x$patient)
  CHECK_EQ(donor.id + (patient.id-1) * length(kBloodCodes), pcode)
  CHECK_EQ(x$desc, sprintf("%s-%s", as.blood.type(x$donor), as.blood.type(x$patient)))
  CHECK_EQ(x$prob, get.blood.code.prob(x$donor) * get.blood.code.prob(x$patient))
  CHECK_EQ(pc.to.pair.type(c(1, 6, 9, 10)), c("S", "S", "O", "R"))
  all.types <- table(pc.to.pair.type(1:16))
  CHECK_EQ(all.types[["R"]], 2)
  CHECK_EQ(all.types[["S"]], 4)
  CHECK_EQ(all.types[["U"]], 5)
  # Check pc.reciprocal (reciprocal pairs)
  rec.id = pc.reciprocal(pcode)
  x2 = subset(kPairs, pc==rec.id)
  CHECK_EQ(c(x2$donor, x2$patient), c(x$patient, x$donor))
}

test.terminology.rpairs <- function() {
  # Checks the sampling functions for pairs
  # Main idea is to match the observed frequency of pairs with specific
  # (pc, pra) , i.e. pair code and PRA sensitivity
  # and test it against the expected value of that frequency.
  CHECK_ERROR(expr=rpairs(0, c(),uniform.pra=T))
  pras = rpra(100, is.uniform=T)
  CHECK_EQ(pras, rep(kUniformPRA, 100))
  pras = rpra(10000, is.uniform=F)
  # Check whether the underlying numbers came from the non-U PRA distribution.
  x = as.vector(table(pras))
  CHECK_GE(chisq.test(x=x, p=kNonUniformPRADistribution)$p.value, 0.01)
  
  marginal.prob <- function(pc, pra, uniform) {
    # Finds the joint probability P(pc, pra, InPool)
    # which is equal to:
    # P(pc, pra) = P(pc) * P(pra) * P(inpool | pc, pra)
    Prob.pc <- kPairs[pc, ]$prob
    Prob.pra <- ifelse(uniform, as.numeric(pra==kUniformPRA), 
                       kNonUniformPRADistribution[which(kNonUniformPRA==pra)])
    Prob.inpool = ifelse(kPairs[pc,]$blood.compatible, pra, 1)
    Prob.pc * Prob.pra * Prob.inpool
  }
  # Calculate the sums.
  sum.marginals <- function(uniform.pra) {
    if(uniform.pra)
      return(sum(marginal.prob(kPairCodes, kUniformPRA, T)))
    return(sum(sapply(kNonUniformPRA, function(i) sum(marginal.prob(kPairCodes, i, F)))))
  }
  
  check.obs.frequency <- function(pairs, arg.pc, arg.pra, uniform.pra) {
    x = subset(pairs, pc==arg.pc & pra==arg.pra)
    ncount = nrow(x)
    prob = marginal.prob(arg.pc, arg.pra, uniform=uniform.pra) / 
          sum.marginals(uniform.pra)
    # Check whether a particular (PC, PRA) has an observed frequency
    # that matches the expected one.  (Uniform Case)
    CHECK_NEAR(ncount / nsamples, prob, tol=3 * sqrt(prob * (1-prob) / nsamples))  
  }
  
  for(uniform.pra in c(T,F)) {
    nsamples = 1000
    all.pairs = rpairs(nsamples, pair.ids=1:nsamples, uniform.pra=uniform.pra)
    random.pc = sample(kPairCodes, 1)
    random.pra = sample(all.pairs$pra, 1)
    check.obs.frequency(all.pairs, random.pc, random.pra, uniform.pra=uniform.pra)
  }
}

test.terminology.edges <- function() {
  # Tests for creating edges from pairs
  # Functions: generate.pairs.edges, map.edges.adjacency
  npairs = 500
  CHECK_ERROR(generate.pairs.edges(pairs=NULL))
  pairs = rpairs(npairs, pair.ids=1:npairs, uniform.pra=T)
  pairs[1,]$pra = 0  # not-sensitive patient (incoming from all compatible pairs)
  pairs[2,]$pra = 1  # very-sensititive patient (no incoming edges)
  pairs[3,]$pra = 0.55
  edges = generate.pairs.edges(pairs)
  compatible.pcs <- function(pc) {
    incoming.pcs <- sapply(kPairCodes, function(i) {
      d2 = kPairs[i,]$donor
      p1 = kPairs[pc,]$patient
      p1 %% d2 == 0
    })
    return(which(incoming.pcs))
  }
  
  # incoming.pcs = (1,5, 15) = PC codes that can give to patient 0
  incoming.pcs = compatible.pcs(pairs[1,]$pc)
  # Check whether all compatible pairs -> patient 0
  incoming.ids <- subset(pairs, pc %in% incoming.pcs)$pair.id
  
  x = subset(edges, pair.id1 %in% incoming.ids & pair.id2 == pairs[1,]$pair.id & self.loop==0,
             select=can.donate)
  # Test #1. Check that non-sensitive patient gets all possible incoming edges
  CHECK_EQ(nrow(x), colSums(x))
  # Check that the sensitive patient has no incoming edges
  x = subset(edges, pair.id2 == pairs[2,]$pair.id & self.loop==0,
             select=can.donate)
  # Test #2. Check that sensitive guy gets no edges (from anyone)
  CHECK_EQ(colSums(x), 0)
  
  # Test #3. Check marginal frequency of edges.
  # Pick a random pair
  random.pair = sample(npairs, 1)
  pra = pairs[random.pair,]$pra
  pc = pairs[random.pair, ]$pc
  id = pairs[random.pair, ]$pair.id
  incoming.pcs = compatible.pcs(pc)
  incoming.ids <- subset(pairs, pc %in% incoming.pcs)$pair.id
  x = subset(edges, pair.id1 %in% incoming.ids & pair.id2 == id & self.loop==0)
  CHECK_MU0(as.vector(x$can.donate), mu0=1-pra)
  
  ## Test map.edges.adjacency
  A = map.edges.adjacency(edges, all.pair.ids=1:npairs)
  CHECK_EQ(A[, 2], rep(0, npairs))
  incoming.pcs = compatible.pcs(pairs[1,]$pc)
  # Check whether all compatible pairs -> patient 0
  incoming.ids <- subset(pairs, pc %in% incoming.pcs)$pair.id
  m = rep(0, npairs)
  m[incoming.ids] <- 1
  CHECK_EQ(A[, 1], m)
  for(i in 1:10) {
    row <- sample(1:nrow(A), 1)
    col <- sample(1:ncol(A), 1)
    p1 <- rownames(A)[row]
    p2 = colnames(A)[col]
    CHECK_EQ(A[row, col], subset(edges, pair.id1==p1 & pair.id2==p2)$can.donate)
  }
}


####    TESTS for rke.R  #####
get.RKE0 <- function() {
  # Creates a testing-rke
  # See a picture at "tests/RKE0.png"
  warning("We can create more complicated graphs here.")
  pairs = subset(kPairs, desc %in% c("A-B", "O-A", "A-A", "B-A", "O-B", "B-B"))
  pairs$pra <- rep(0, nrow(pairs))
  pairs$pair.id <- 1:nrow(pairs)
  pairs$hospital <- rep(1, nrow(pairs))
  rke = empty.rke()
  rke$pairs = pairs
  rke <- rke.update.new.pairs(rke, keep.edges=c())
  edg <- data.frame(from=1:6, to=c(6, 3, 4, 2, 1, 5))
  give.to <- function(i) {
    x = c(subset(edg, from==i)$to)
    if(i == 3) return(c(x, 5))
    if(i == 5) return(c(x, 3))
    return(x)
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
  rke = get.RKE0()
  cycles2 = rke.cycles(rke, include.3way=F)
  cycles3 = rke.cycles(rke, include.3way=T)
  CHECK_EQ(nrow(cycles2), 1)
  CHECK_EQ(nrow(cycles3), 3)
}

#     Tests for MATCHING.
test.matching <- function() {
  rke <- get.RKE0()
  m = max.matching(rke, include.3way=F)
  CHECK_EQ(m$utility, 2)
  m = max.matching(rke, include.3way=T)
  CHECK_EQ(m$utility, 6)
}

####    Testing for mechanisms.R
test.rCM <- function() {
  source("mechanisms.R")
  rke.pool <- rrke.pool(m=3, n=30, uniform.pra=T)
  kpd <- kpd.create(rke.pool, "ttt")
  U = Run.Mechanism(kpd, "rCM")
  m = max.matching(rke.pool$rke.all)
  U2 = get.hospitals.utility(rke.all=rke.pool$rke.all, matched.ids=m$match$pair.id)
  # All truthful matching
  CHECK_EQ(U, U2)
}





