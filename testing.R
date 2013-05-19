## Unit tests.
TEST.SETS.DISJOINT = function(x, y, str="n/a") {
  if(length(intersect(x,y))>0)
    stop(sprintf("[TEST FAIL]...Sets x,y not disjoint : %s", str))
  return(T)
}
TEST.SETS.EQUAL = function(x, y, str="n/a") {
  throw = function() {
    stop(sprintf("[TEST FAIL]...Sets are not equal : %s", str))
  }
  if(length(x) != length(y))
    throw()
  if(length(setdiff(x,y))>0)
    throw()
  return(T)
}
TEST.LISTS.EQ = function(x,y, tol=0, str="n/a") {
  throw = function() {
    stop(sprintf("[TEST FAIL]...Lists x,y not equal at this tol level : %s", str))
  }
  s = sum(which(abs(x-y)>tol))
  if(s >0) throw()
  return(T)
}
## Test whether two RKE's are the same (NOT exact)
## TEST.RKE.EQ(rrke(10), rrke(10))   should fail
TEST.RKE.EQ = function(x,y, str="n/a") {
  throw = function(astr="n/a") {
    stop(sprintf("[TEST FAIL]..RKEs are not equal : %s : %s", astr, str))
  }
  if(get.size(x) != get.size(y))
    throw("size")
  
  library(igraph)
  gx = graph.adjacency(x$P * x$B, mode="undirected")
  gy = graph.adjacency(y$P * y$B, mode="undirected")
  
  if(vcount(gx)!= vcount(gy)) throw("vcount")
  if(ecount(gx) != ecount(gy)) throw("ecount")
  if(length(get.diameter(gx)) != length(get.diameter(gy))) throw("diameter")
  if(mean(alpha.centrality(gx)) != mean(alpha.centrality(gy))) throw("centrality")
  if(mean(degree.distribution(gx)) != mean(degree.distribution(gy))) throw("degree")
  
  return(T)
  
}
test.repeat = function(test, args, trials) {
  for(t in 1:trials) {
    res = test(args)
    print(sprintf("t=%d/%d:   [OK]", t, trials))
  }
  
}

### Tests compute.ir.constraints from mechanisms.R
test.ir.constraints = function(args) {
  m = args$m
  n= args$n
  
  rke.list = rrke.many(m=m, n=n)
  ir = compute.ir.constraints(rke.list, types=c("S","R"))
  to.vector = function(found.pcs) {
    x = rep(0, length(Pair.Codes))
    for (i in unique(found.pcs)) {
      no.pc = length(which(found.pcs==i))
      x[i]=no.pc
    }
    return(x)
  }
  ##  S-test
  for(hid in 1:length(rke.list)) {
    rke = rke.list[[hid]]
    rs = get.subgraph(rke, "S")
    m = max.matching(rs)
    pcs = to.vector(  rs$pc[m$matching$matched.ids]   )
    
    if(! TEST.LISTS.EQ(pcs, ir$S[[hid]]) ) return(F)
  }  
  ##  R-test
  for(hid in 1:length(rke.list)) {
    rke = rke.list[[hid]]
    rs = get.subgraph(rke, "R")
    m = max.matching(rs)
    pcs = to.vector(  rs$pc[m$matching$matched.ids]  )
    if(! TEST.LISTS.EQ(pcs, ir$R[[hid]])) return(F)
  }   
  return(T)
}
test.init.mechanism = function(args) {
  m = args$m
  n = args$n 
  rke.list = rrke.many(m=m, n=n)
# part 1. Test truthful strategies
  x = init.mechanism(rke.list, paste(rep("t", m), collapse=""))
  for(hid in 1:m) {
    rke = rke.list[[hid]]
    ## Check if the utility
    TEST.LISTS.EQ(0, x$util[hid,],str=" t - utility test")
    TEST.RKE.EQ(rke, x$rke.list[[hid]], str="t-rke test")
  }
  ## part 2. test when all canonically deviate
  x = init.mechanism(rke.list,paste(rep("c", m), collapse=""))
  for(hid in 1:m) {
    rke = rke.list[[hid]]
    ## Check if the utility
    m = max.matching(rke)
    TEST.LISTS.EQ(m$matching$utility, x$util[hid,], str="c- utility test")
    TEST.RKE.EQ(x$rke.list[[hid]], remove.pairs(rke, m$matching$matched.ids),
                str="c- rke test")
  }
  return(T)
  
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
