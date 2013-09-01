# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
source("../r-toolkit/checks.R")

test.repeat = function(test, args, trials) {
  for(t in 1:trials) {
    res = test(args)
    print(sprintf("t=%d/%d:   [OK]", t, trials))
  }  
}

test.terminology <- function() {
  x = rpra(10000, is.uniform=F)
  test.out = chisq.test(table(x), p=kNonUniformPRADistribution)
  CHECK_TRUE(test.out$p.value > 0.01,
             msg=c(sprintf("df=%d", test.out$parameter), test.out$statistic))
  
  rke = rrke(1000, uniform.pra=T)
  CHECK_rke(rke)
  estim.var = kUniformPRA * (1-kUniformPRA) / 1000
  CHECK_NEAR(x=mean(rke$edges$pra.compatible), y=1-kUniformPRA,
             tol= 2 * sqrt(estim.var), msg="Mean obs. PRA")
  CHECK_SETEQ(subset(rke$pairs, pair.type=="U")$blood.compatible, c(0),
              msg="U pairs should not be blood-type compatible")
}




