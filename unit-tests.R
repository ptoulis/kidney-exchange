## Unit tests.
# 1.  Sampling RKE. Distribution/expected values of ABO groups?
# 2.  Inconsistencies in RKE graphs (e.g. self-loops, matches correct, missing matches?)
# 3.  Shuffling in Gurobi max matching?
# 4. Regular matching?
# 5. Removing pairs from RKE objects.
#  6. Is pooling of RKE's deterministic ?  Returns the same RKE 
## 7.  if PRA=0  then what results you expect?
test.pooling <- function() {
    rkes = rrke.many(3, n=20)
    r1 = pool.rke(rkes)
    for(i in 1:10) {
        r2 = pool.rke(rkes)
        match = (get.size(r1)==get.size(r2)) &&
            (nrow(get.edgelist(r1$graph)) == nrow(get.edgelist(r2$graph)));
        if(!match) {
            save(rkes, file="debug/pooling-problem.Rdata")
            stop("Test FAIL.")
        }
    }
    print("Test PASS")
}

perf.test <- function() {
    print("Running Performance Tests")
    print("1. Time to sample 300-size and 600-size graphs:")
    print(system.time(rrke(300)))
    print(system.time(rrke(600)))
    print("2. Experiment 3 profile (3 hospitals, 50 sims, size=20)")
    prof.file= "debug/perftest.prof"
    Rprof(filename=prof.file)
    E3(ns=1)
    Rprof(NULL)
    summaryRprof(filename=prof.file)
}


test.pairs = function() {
    rke = rrke(100)
    A = get.model.A(rke)
    a0 =  filter.edges.by.type(rke, A, "O", "U")
    a1 = filter.edges.by.type(rke, A, "O", "O")
    a2 = filter.edges.by.type(rke, A, "O", "R")
    a3 = filter.edges.by.type(rke, A, "O", "S")
    a4 = filter.edges.by.type(rke, A, "R", "R")
    a5 = filter.edges.by.type(rke, A, "S", "S")
    all = c(a0,a1,a2,a3,a4,a5)
    print(ncol(A))
    print(length(all))
    original = 1:ncol(A)
    d= setdiff(original, all )
    print(sprintf("Edge lost %d", d[1]))
    ids = which(A[,d[1]]==1)
    print(sprintf("Difference %d ", length(d)))
    print(ids)
    pcs = rke$pc[ids]
    print(pcs)
    print(pair.code.to.pair(pcs))
}


test.xCM = function() {
  source("rke2.R")
  source("mechanisms.R")
  ##  Now test mechanism. First test if rCM >= xCM
  trials = 100
  for(i in 1:trials) {
    rke.list = rrke.many(m=3, n=20)
    sumr = sum(rCM(rke.list, "ttt"))
    sumx = sum(xCM(rke.list, "ttt"))
    if(sumr >= sumx) {
      print("[PASS]...")
      
    } else {
      print("[FAIL].. Saving state.")
      save(rke.list, file="debug/xcm-fail.Rdata")
      break()
    }
  }
}
