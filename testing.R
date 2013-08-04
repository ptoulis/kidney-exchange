## Unit tests.
stop.now <- function(x, y, msg) {
  print(sprintf("[FAIL]: %s", paste(msg, collapse=" : ")))
  str(x, width=80)
  str(y, width=80)
  stop("")
}

CHECK_TRUE <- function(x, msg="n/a") {
  if (!x)
    stop.now(x, "n/a (ignore)", c("Logical test failed", msg))
}

CHECK_EQ <- function(x, y, msg="n/a") {
  if (!all(x==y)) stop.now(x, y, c("Arrays not equal.", msg))
  return (TRUE)  
}

CHECK_GE <- function(x, y, msg="n/a") {
  if (any(x < y))
    stop.now(x, y, c("Some values are > y", msg))
}

CHECK_SETEQ <- function(x, y, msg="n/a") {
  if (!setequal(x, y))
    stop.now(x,y, c("Sets not equal", msg))
}

CHECK_DISJOINT <- function(x, y, msg="n/a") {
  if(length(intersect(x, y)) > 0)
    stop.now(x, y, c("Sets are not disjoint", msg))
}

CHECK_UNIQUE <- function(x, msg="n/a") {
  if (any(duplicated(x)))
    stop.now(x, "na", c("Duplicates exist", msg))
}

CHECK_MEMBER <- function(element, y, msg="n/a") {
  if (!all(is.element(element, y))) {
    stop.now(element, y, sprintf("Element not member of list : %s", msg))
  }
  return(TRUE)
}

CHECK_NEAR <- function(x, y, tol=1e-2, msg="n/a") {
  if (any(abs(x-y) > tol))
    stop.now(x, y, c(sprintf("Failed at tolerance %.6f", tol), msg))
}

CHECK_INTERVAL <- function(x, min, max, msg="n/a") {
  if(any(x < min) | any(x > max))
    stop.now(x, c(min, max), c("Not in interval", msg))
}

CHECK_EXCEPTION <- function(expr, msg="n/a") {
  out <- tryCatch(eval(expr), error=function(err) { return("ERROR") })
  CHECK_TRUE(out == "ERROR", msg=c("Did not throw exception", msg))
}
bootstrap.mean = function(x) sd(replicate(1000, { mean(sample(x, replace=T)) }))
CHECK_MU0 <- function(x, mu0, msg="n/a") {
  # Will bootstrap x to see if mean(x) = mu0
  boot.se = bootstrap.mean(x)
  CHECK_INTERVAL(mu0, min=mean(x) - 2 * boot.se, max=mean(x) + 2*boot.se,
                 msg=c("Not in bootstrap inteval", msg))
}

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




