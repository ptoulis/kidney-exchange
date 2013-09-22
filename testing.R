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




