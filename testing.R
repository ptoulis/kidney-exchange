# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
rm(list=ls())
# Loads everything.
source("experiments.R")

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

test.matching <- function() {
  ntriangles = rpois(1, lambda=20)
  rke = get.RKE0(ntriangles)
  loginfo(sprintf("Testing with %d triangles..", ntriangles))
  m = max.matching(rke, include.3way=F)
  CHECK_EQ(m$utility, ifelse(ntriangles %% 2 ==0, ntriangles, ntriangles-1))
  m = max.matching(rke, include.3way=T)
  CHECK_EQ(m$utility, 3 * ntriangles)
}

test.matching.randomness <- function(nsamples=10) {
  nhospitals = 2
  rke.all <- rrke.pool(m=nhospitals, n=40, uniform.pra=T)$rke.all
  rke.all$A <- 1 * (rke.all$A + t(rke.all$A) > 0)
  # create symmetric hospitals
  U = matrix(0, nrow=0, ncol=nhospitals)
  Unorandom <- matrix(0, nrow=0, ncol=nhospitals)
  pb = txtProgressBar(style=3)
  for (i in 1:nsamples) {
    m = max.matching(rke.all, include.3way=F)
    m2 <- max.matching(rke.all, include.3way=F, randomize.matching=F)
    
    U <- rbind(U, get.matching.hospital.utilities(m, nhospitals))
    Unorandom <- rbind(Unorandom,  get.matching.hospital.utilities(m2, nhospitals))
    print(sprintf("%s", paste(round(colMeans(U), 3), collapse=", "))) 
    print(sprintf("%s", paste(round(colMeans(Unorandom), 3), collapse=", "))) 
    # setTxtProgressBar(pb, value=i/nsamples)
  }
  print("Checking all maximum matchings.")
  CHECK_SETEQ(rowSums(U), c(sum(U[1, ])))
  print("Observed counts.")
  print(U)
  print("Colsums")
  print(colSums(U))
  print("Test")
  print(chisq.test(colSums(U)))
}



#     Tests for MATCHING.
check.rke.strategy <- function(rke, strategy.out, strategy.str) {
  rke.hide <- rke.keep.pairs(rke, strategy.out$hide)
  rke.report <- rke.keep.pairs(rke, strategy.out$report)
  # union of hide + report gives the original one.
  CHECK_TRUE(rke.equal(rke.add(rke.hide, rke.report), rke))
  CHECK_DISJOINT(strategy.out$hide, strategy.out$report)
  CHECK_SETEQ(union(strategy.out$hide, strategy.out$report), rke.pair.ids(rke))
  # strategy-specific considerations
  if (strategy.str=="t") {
    # no hiding
    CHECK_EQ(length(strategy.out$hide), 0)
    CHECK_TRUE(rke.equal(rke$report, rke))
  } else if (strategy.str=="c") {
    m = max.matching(rke)
    # matched ids should be same in # with a max matching
    CHECK_EQ(length(get.matching.ids(m)), length(strategy.out$hide))
  } else {
    nAB <- nrow(subset(rke$pairs, desc=="A-B"))
    nBA <- nrow(subset(rke$pairs, desc=="B-A"))
    CHECK_EQ(length(strategy.out$hide), min(nAB, nBA))
  }
}

test.play.strategy <- function() {
  rke <- rrke(20)
  CHECK_ERROR(play.strategy(rke, strategy.str="p", include.3way=F),
              msg="Not valid strategy")
  for (s in c("t", "c", "r")) {
    out <- play.strategy(rke, s)
    check.rke.strategy(rke, out, s)
  }
}

test.play.strategies <- function() {
  pool <- rrke.pool(m=3, n=20, uniform.pra=T)
  # wrong no. of strategies.
  CHECK_ERROR(play.strategies(rke.list=pool$rke.list, strategy.str="tttttt"))
  
  strategies <- sample(c("t", "c", "r"), size=3, replace=T)
  strategy.str = paste(strategies, collapse="")
  out = play.strategies(rke.list=pool$rke.list, strategy.str=strategy.str)
  hids = rke.list.hospital.ids(rke.list=pool$rke.list)
  for(i in hids) {
    loginfo(sprintf("Checking hospital %d strategy=%s", i, strategies[i]))
    check.rke.strategy(pool$rke.list[[i]], out[[i]], strategies[i])
  }
}

test.rrke.pool <- function(nsamples=10) {
  # Check whether pooling through rrke.pool()
  # and using rrke() yield the same observed frequencies for pairs.
  rrke.table <- c()
  pool.table <- c()
  kpd.table <- c()
  nsize = 30
  pb <- txtProgressBar(style=3)
  for(i in 1:nsamples) {
    pool <- rrke.pool(m=3, n=nsize / 3, uniform.pra=T)
    kpd <- kpd.create(pool, strategy.str="ttt")
    rke <- rrke(nsize, uniform.pra=T)
    if(length(pool.table)==0) {
      pool.table = table(pool$rke.all$pairs$desc)
    } else {
      pool.table = pool.table + table(pool$rke.all$pairs$desc)
    }
    
    if (length(rrke.table)==0) {
      rrke.table = table(rke$pairs$desc)
    } else {
      rrke.table = rrke.table + table(rke$pairs$desc)
    }
    if(length(kpd.table) == 0) {
      kpd.table <- table(kpd$real.pool$rke.all$pairs$desc)
    } else {
      kpd.table <- kpd.table+ table(kpd$real.pool$rke.all$pairs$desc)
    }
      
    setTxtProgressBar(pb, value=i/nsamples)
  }
  rmu <- 100 * rrke.table / (nsize * nsamples)
  pmu <- 100 * pool.table / (nsize * nsamples)
  kmu <- 100 * kpd.table  / (nsize  * nsamples)
  CHECK_NEAR(kmu, pmu, tol=1e-3)
  print("Comparing observed frequencies")
  print(round(rmu, 4))
  print(round(pmu, 4))
  print(round(kmu, 4))
  print(chisq.test(rmu, pmu))
  print(chisq.test(rmu, kmu))
}

test.kpd.create <- function() {
  pool = rrke.pool(3, 20, T)
  strategies <- sample(c("t", "c", "r"), size=3, replace=T)
  strategy.str = paste(strategies, collapse="")
  loginfo(sprintf("Strategies= %s", strategy.str))
  kpd <- kpd.create(rke.pool=pool, strategy.str=strategy.str)
  for(hid in 1:length(kpd$reported.pool$rke.list)) {
    rke.report <- kpd$reported.pool$rke.list[[hid]]
    rke.real <- kpd$real.pool$rke.list[[hid]]
    CHECK_TRUE(rke.equal(rke.real, pool$rke.list[[hid]]))
    out = list(report=rke.pair.ids(rke.report),
               hide=setdiff(rke.pair.ids(rke.real), rke.pair.ids(rke.report)))
    check.rke.strategy(rke.real, out, strategies[hid])
  }
}


####    Testing for mechanisms.R
test.rCM <- function() {
  source("mechanisms.R")
  nhospitals = 3
  rke.pool <- rrke.pool(m=nhospitals, n=30, uniform.pra=T)
  kpd <- kpd.create(rke.pool, "ttt")
  U1 = get.matching.hospital.utilities(Run.Mechanism(kpd, "rCM", include.3way=F), nhospitals)
  match = max.matching(rke.pool$rke.all)
  U2 = get.matching.hospital.utilities(match, nhospitals)
  
  # All truthful matching
  CHECK_EQ(U1, U2)
}



create.utility.matrix <- function(m) {
  M <- matrix(0, nrow=3, ncol=3)
  rownames(M) <- c("H1", "H2", "H3")
  colnames(M) <- c("Mech", "Internal", "Total")
  Um <- get.matching.hospital.utilities(m$mech.matching, 3)
  Ut <- get.matching.hospital.utilities(m$total.matching, 3)
  for(i in 1:3) {
    M[i, 1] = Um[i]
    M[i, 2] = get.matching.utility(m$internal.matchings[[i]])
    M[i, 3] = Ut[i]
  }
  return(M)
}




create.Ronly.rke <- function(nAB, nBA, hid, start.pair.id, uniform.pra) {
  all.pairs <- subset(kPairs, desc=="none")
  for (i in 1:nAB)
    all.pairs <- rbind(all.pairs, subset(kPairs, desc=="A-B"))
  for(i in 1:nBA) 
    all.pairs <- rbind(all.pairs, subset(kPairs, desc=="B-A"))
  
  N <- nrow(all.pairs)
  all.pairs$hospital <- rep(hid, N)
  all.pairs$pair.id <- start.pair.id + seq(0, N-1)
  all.pairs$pra <- rpra(N, is.uniform=uniform.pra)
  
  rke = empty.rke()
  rke$pairs <- all.pairs
  rke <- rke.update.new.pairs(rke, keep.edges=NULL)
  return(rke)
}

rbipartite.allocation <- function(a, b, y) {
  ## Consider the following problem:
  #   AB     BA
  #   (a)    (b+y)
  #
  # i.e. need to allocated a AB pairs to b+y BA pairs.
  # How #many  BA pairs from the y in total will be used?
  # Define this as N(b, y; a).  We know E N =  a * y / (y + b)
  # What about the distribution?
  if (y == 0) return(0)
  if (a==0) return(0)
  # know that a, y != 0
  if(b==0) return(1)
  
  X <- rbinom(1, size=1, prob = y / (y+b))
  if(X == 1) return(rbipartite.allocation(a-1, b, y-1) + 1)
  return(rbipartite.allocation(a-1, b-1, y))
}

Nmatchings <- function(a, b, x, y) {
  # Assume x AB, y BA pairs    and  x AB, y BA for some Hospital Hi
  # Returns matches for Hi assuming random bipartite matching.
  CHECK_GE(b + y, a + x)
  N <- x + replicate(5000, { rbipartite.allocation(a + x, b, y) })
  return(N)
}

CHECK_Rpairs <- function(Rpairs) {
  CHECK_SETEQ(names(Rpairs), c("AB", "BA"))
  CHECK_EQ(nrow(Rpairs), 3)
}

mech.theor.matchings.Rdev <- function(mech, h3.strategy,
                                      Rpairs) {
  CHECK_MEMBER(h3.strategy, c("t", "r"))
  CHECK_Rpairs(Rpairs)
  a = sum(head(Rpairs$AB, 2))
  b = sum(head(Rpairs$BA, 2))
  x = Rpairs$AB[3]
  y = Rpairs$BA[3]
  
  if(mech=="rCM") {
    if(h3.strategy == "t") {
      return(Nmatchings(a, b, x, y))
    } else {
      ## strategy = r
      N = Nmatchings(a, b, x=0, y)
      y.remainder = y - N
      internal <-2 * sapply(y.remainder, function(i) min(i, x))
      return(N + internal)
    }
  } else if (mech=="Bonus") {
    short.side = head(apply(Rpairs, 1, min), 2)
    pairtype = head(apply(Rpairs, 1, which.min), 2)
    s = sum(short.side[pairtype==1])
    CHECK_EQ(s, sum(short.side[pairtype==2]))
    a = a - 2 * s
    b = b - 2 * s
    print(sprintf("Bonus mechanism matches 4 * %d = %d internally", s , 4 * s))
    if(h3.strategy == "t") {
      return(Nmatchings(a, b, x, y))
    } else {
      ## strategy = r
      N = Nmatchings(a, b, x=0, y)
      y.remainder = y - N
      internal <-2 * sapply(y.remainder, function(i) min(i, x))
      return(N + internal)
    }
  }
}

example.Rpairs <- function() {
  return(data.frame(AB=c(40, 10, 15), 
                    BA=c(10, 40, 30)))
}

test.Rdeviation <- function(ntrials=10, Rpairs=example.Rpairs(),
                            mech="rCM", pra.value=0) {
  # Tests only on the R subgraph, whether the R-deviation works.
  # Creates a R-only graph.
  #
  #  Example, nAB=30   nBA=10
  #  H1=(30, 10)   H2=(10, 30)  H3=(h3.nAB, h3.nBA)
  CHECK_Rpairs(Rpairs)
  H3.truth.util <- c()
  H3.Rdev.util <- c()
  kUniformPRA <<- pra.value
  
  print("Computing theoretical differences.")
  ## Compute theoretical 
  theoretical.true.matches <-  mech.theor.matchings.Rdev(mech=mech, h3.strategy="t",
                                                         Rpairs=Rpairs)
  theoretical.rdev.matches <- mech.theor.matchings.Rdev(mech=mech, h3.strategy="r",
                                                        Rpairs=Rpairs)
  for(trials in 1:ntrials) {
    print(sprintf("Running mechanism %s", mech))
    rke.pool <- list(rke.list=list(), rke.all=empty.rke())
    
    r1 <- create.Ronly.rke(Rpairs$AB[1], Rpairs$BA[1], hid=1, 1, uniform.pra=T)
    r2 <- create.Ronly.rke(Rpairs$AB[2], Rpairs$BA[2], hid=2, 100, uniform.pra=T)
    r3 <- create.Ronly.rke(Rpairs$AB[3], Rpairs$BA[3], hid=3,  200, uniform.pra=T)
    
    print("Sampled PRA")
    print(r1$pairs$pra)
    # print(r2$pairs$pra)
    rke.pool$rke.list[[1]] <- r1
    rke.pool$rke.list[[2]] <- r2
    rke.pool$rke.list[[3]] <- r3
    
    for(i in 1:3) {
      rke = rke.pool$rke.list[[i]]
      m = max.matching(rke)
      print(sprintf("H-%d: AB=%d BA=%d : Internally matching %d pairs",
                    i, nrow(subset(rke$pairs, desc=="A-B")), nrow(subset(rke$pairs, desc=="B-A")),
                    get.matching.utility(m)))
      rke.pool$rke.all <- rke.add(rke.pool$rke.all, rke.pool$rke.list[[i]])
    }
    
    kpd.t <- kpd.create(rke.pool, "ttt")
    kpd.r <- kpd.create(rke.pool, "ttr")
    out.t = Run.Mechanism(kpd.t, mech=mech, include.3way=F)
    out.r = Run.Mechanism(kpd.r, mech=mech, include.3way=F, verbose=F)
    
    Utotal.t <- get.matching.hospital.utilities(out.t$total.matching, 3)
    Utotal.r <- get.matching.hospital.utilities(out.r$total.matching, 3)
    print(sprintf("ttt: %s | total = %d", paste(Utotal.t, collapse=", "), sum(Utotal.t)))
    print(create.utility.matrix(out.t))
    
    print(sprintf("ttr: %s | total = %d", paste(Utotal.r, collapse=", "), sum(Utotal.r)))
    print(create.utility.matrix(out.r))
    
    H3.Rdev.util <- c(H3.Rdev.util, Utotal.r[3])
    H3.truth.util <- c(H3.truth.util, Utotal.t[3])
    print.stats <- function(vec, str) {
      mu <- mean(vec)
      se <- bootstrap.mean(vec)
      print(sprintf("%s : %.2f, CI=[%.2f, %.2f]", str, mu, mu-2*se, mu+2*se))
    }
    print.stats(H3.Rdev.util, "R-deviation estimate")
    print.stats(theoretical.rdev.matches, str="R-dev theoretical: ")
    print("******")
    print.stats(H3.truth.util, "Truthful estimate")
    print.stats(theoretical.true.matches, str="Truth theoretical: ")
    print("*****")
    print.stats(H3.Rdev.util - H3.truth.util, "Empirical Difference")
    print.stats(theoretical.rdev.matches - theoretical.true.matches, "Theoretical difference")
  }
}


