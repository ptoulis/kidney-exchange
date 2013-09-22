# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
source("../r-toolkit/checks.R")

test.repeat = function(test, trials) {
  for(t in 1:trials) {
    res = test()
    print(sprintf("t=%3d/%3d:   [OK]", t, trials))
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

test.rpairs <- function() {
  # Checks the sampling functions for pairs.
  pras = rpra(100, is.uniform=T)
  CHECK_EQ(pras, rep(kUniformPRA, 100))
  pras = rpra(10000, is.uniform=F)
  # Check whether the underlying numbers came from the non-U PRA distribution.
  x = as.vector(table(pras))
  CHECK_GE(chisq.test(x=x, p=kNonUniformPRADistribution)$p.value, 0.01)
}




