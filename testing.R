## Unit tests.
source("lib.R")
equal.sets = function(x,y) {
  if(length(x) != length(y)) return(F)
  return(length(setdiff(x,y))==0)
}
is.subset = function(bigger, smaller) {
  xandy = intersect(bigger, smaller)
  return(equal.sets(xandy, smaller))
}

TEST.SETS.DISJOINT = function(x, y, str="n/a") {
  if(length(intersect(x,y))>0)
    stop(sprintf("[TEST FAIL]...Sets x,y not disjoint : %s", str))
  return(T)
}
## Testing set-equality
TEST.SETS.EQ = function(x, y, str="n/a") {
  x = unique(x)
  y = unique(y)
  throw = function() {
    stop(sprintf("[TEST FAIL]...Sets are not equal : %s", str))
  }
  if(length(x) != length(y))
    throw()
  if(length(setdiff(x,y))>0)
    throw()
  return(T)
}
TEST.LISTS.EQ = function(x,y, str="n/a", tol=0) {
  throw = function() {
    stop(sprintf("[TEST FAIL]...Lists x,y not equal at this tol level : %s", str))
  }
  s = sum(which(abs(x-y)>tol))
  if(s >0) throw()
  return(T)
}
TEST.LISTS.GEQ = function(x,y, str="n/a") {
  throw = function() {
    stop(sprintf("[TEST FAIL]...Does not hold x >= y : %s", str))
  }
  s = sum(which(x-y<0))
  if(s >0) throw()
  return(T)
}
TEST.SUBSET = function(smaller, bigger, str="n/a)") {
  smaller = unique(smaller)
  bigger = unique(bigger)
  throw = function() {
    stop(sprintf("[TEST FAIL]...x is not subset of y : %s", str))
  }
  if(! is.subset(bigger=bigger, smaller=smaller))
    throw()
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
  gx = rke.to.igraph(x)
  gy = rke.to.igraph(y)
  
  if(vcount(gx)!= vcount(gy)) throw("vcount")
  if(ecount(gx) != ecount(gy)) throw("ecount")
  if(length(get.diameter(gx)) != length(get.diameter(gy))) throw("diameter")
  if(mean(alpha.centrality(gx)) != mean(alpha.centrality(gy))) throw("centrality")
  if(mean(degree.distribution(gx)) != mean(degree.distribution(gy))) throw("degree")
  
  return(T)
  
}
TEST.MEAN.EQ = function(x, mu0, str="n/a") {
  sd = bootstrap.mean(x)
  mu = mean(x)
  ci = c(mu - 2 *sd, mu + 2*sd)

  if(mu0 < ci[1] || mu0 > ci[2])
    stop(sprintf("[TEST FAIL] Hypothesis test mu=mu0  mu0=%.3f  CI=[%.3f, %.3f]: %s", 
                 mu0, ci[1], ci[2], str))
  if(sd> mu/2)
    warning("SE is probably too high. Try increasing the sample size?")
  return(T)
}
TEST.BOOL = function(x, str="n/a") 
{
  if(! x ) 
    stop(paste("[TEST FAIL]...test.bool: %s", str) )
  return(T)
}



test.repeat = function(test, args, trials) {
  for(t in 1:trials) {
    res = test(args)
    print(sprintf("t=%d/%d:   [OK]", t, trials))
  }
  
}




## Tests for lib.R
test.rpra = function(args) {
  throw=function(str) stop(sprintf("[TEST FAIL]...rpra() hypothesis test fail : %s", str))
  
  n = args$n
  ## sample PRA uniformly
  for(pra.mode in c(T,F)) {
    pra1 = rpra(n,is.uniform=pra.mode)
    pra2 = pra1
    
    Q = rpra.matrix(pra1, pra2);
    if(nrow(Q)!=ncol(Q)) throw("square matrix")
    if(nrow(Q)!=n) throw("correct size")
    ## Pick some ij  pair
    pair = sample(1:n, 2, replace=F)
    i = pair[1]
    j = pair[2]
    reps =replicate(2000, { P = rpra.matrix(pra1, pra2); P[i,j]})
    mu = mean(reps)
    se = bootstrap.mean(reps)
    ci=c(mu-2*se, mu+2*se)
    pci = pra1[i]
    pcj = pra2[j]
    
    p.theor = (1-pci) * (1-pcj)
    if(p.theor< ci[1] || p.theor>ci[2])
      throw(sprintf("Uniform? %s Pij (expected)=%.3f  found=[%.3f, %.3f]",pra.mode,
                    p.theor,  ci[1], ci[2]))
    
  }
  return(T)
}

## Tests for  rke.R
test.rrke = function(args) {
  throw=function(str) stop(sprintf("[TEST FAIL]...rrke() : %s", str))
  
  n = args$n
  trials = 100
  reps= replicate(trials, {
    rke = rrke(n)
    
    oa = length( filter.pairs.by.donor.patient(rke,"O", "A") )
    ao = length( filter.pairs.by.donor.patient(rke,"A", "O") )
    return(ifelse(ao > 0,  oa/ao, NA))
  });
  
  Pc = rpra(1)
  if(sum(is.na(reps)) > 0.2 * trials) 
    stop("Too many NA's. Please try to increase n = size of RKE")
  
  reps = reps[! is.na(reps)]
  # 1. test whether UD pairs are 1/pc  more frequent than the reciprocal OD pairs.
  TEST.MEAN.EQ(reps, mu0=Pc, str="rrke() -- testing OD/UD=pc")
  # 2. dumb test. Check whether the size is ok.
  TEST.LISTS.EQ(get.size(rrke(n)), n,str="Total size")
  
  
  
}
test.pool.rke = function(args) {
  throw=function(str) stop(sprintf("[TEST FAIL]...pool.rke() : %s", str))
  rke.list= rrke.many(m=args$m, n=args$n, uniform.pra=T)
  
  rke.all = pool.rke(rke.list)
  if(length(unique(rke.all$hospital)) != length(rke.list))
    throw("hospital #'s are different!")
  
  for(i in 1:length(rke.list)) {
    hpairs = length(which(rke.all$hospital==i))
    if(hpairs!= get.size(rke.list[[i]]))
      throw("Hospital i does not have the same # of pairs.")
    
    rke.h.pairs=  get.hospital.pairs(rke.all, hid=i)
    rke.h.edges = length(get.internal.edges(rke.all, rke.h.pairs))
    if(rke.h.edges != length(rke.edges(rke.list[[i]])))
      throw("Hospital i does not have the same # of edges")
  }
  
  return(T)
  
}
## Test all filter.*.by.* functions
test.filters = function(args) {
  load(file="tests/rke-20.Rdata")
  #  rke  is a 20-pair hospital
  aa = 1
  oo = 2
  ab = 2
  ba = 2
  ob = 0
  ao = 8
  
  TEST.LISTS.EQ(aa, length(filter.pairs.by.donor.patient(rke,"A", "A") ), str="filter.pairs.by.dp" )
  TEST.LISTS.EQ(oo, length(filter.pairs.by.donor.patient(rke,"O", "O") ), str="filter.pairs.by.dp"  )
  TEST.LISTS.EQ(ab, length(filter.pairs.by.donor.patient(rke,"A", "B") ), str="filter.pairs.by.dp"  )
  TEST.LISTS.EQ(ba, length(filter.pairs.by.donor.patient(rke,"B", "A") ), str="filter.pairs.by.dp"  )
  TEST.LISTS.EQ(ob, length(filter.pairs.by.donor.patient(rke,"O", "B") ), str="filter.pairs.by.dp"  )
  TEST.LISTS.EQ(ao, length(filter.pairs.by.donor.patient(rke,"A", "O") ), str="filter.pairs.by.dp" )
  
  ud = 12
  r = 4
  od = 1
  s = 3
  TEST.LISTS.EQ(ud, length(filter.pairs.by.type(rke, "U")), str="filter.pairs.by.type") 
  TEST.LISTS.EQ(od, length(filter.pairs.by.type(rke, "O")), str="filter.pairs.by.type") 
  TEST.LISTS.EQ(r, length(filter.pairs.by.type(rke, "R")), str="filter.pairs.by.type") 
  TEST.LISTS.EQ(s, length(filter.pairs.by.type(rke, "S")), str="filter.pairs.by.type") 
  
  ob.edges = 0
  oa.edges = 9
  ab.o.edges = 0
  ab.edges = 3
  ba.edges = 3
  TEST.LISTS.EQ(ob.edges, length(filter.edges.by.donor.patient(rke, "O","B")), str="filter.edges.by.dp")
  TEST.LISTS.EQ(oa.edges, length(filter.edges.by.donor.patient(rke, "O","A")), str="filter.edges.by.dp") 
  TEST.LISTS.EQ(ab.o.edges, length(filter.edges.by.donor.patient(rke, "AB","O")), str="filter.edges.by.dp" )
  TEST.LISTS.EQ(ab.edges, length(filter.edges.by.donor.patient(rke, "A","B")), str="filter.edges.by.dp") 
  TEST.LISTS.EQ(ba.edges, length(filter.edges.by.donor.patient(rke, "B","A")), str="filter.edges.by.dp") 
  
  re = 3
  ue = 6
  oe=9
  se = 4
  sr = 0
  os = 3
  TEST.LISTS.EQ(re, length(filter.edges.by.type(rke, "R", "R")), str="filter.edges.by.type") 
  TEST.LISTS.EQ(ue, length(filter.edges.by.type(rke, "U", "*")), str="filter.edges.by.type") 
  TEST.LISTS.EQ(oe, length(filter.edges.by.type(rke, "*", "O")), str="filter.edges.by.type") 
  TEST.LISTS.EQ(se, length(filter.edges.by.type(rke, "S", "*")), str="filter.edges.by.type") 
  TEST.LISTS.EQ(sr, length(filter.edges.by.type(rke, "R", "S")), str="filter.edges.by.type") 
  TEST.LISTS.EQ(os, length(filter.edges.by.type(rke, "S", "O")), str="filter.edges.by.type") 
  
  
}
test.get.model.A = function(args) {
  load(file="tests/rke-20.Rdata")
  A = get.model.A(rke)
  nedges = 14
  # 1. Test no. of edges/pairs
  TEST.LISTS.EQ(ncol(A), 14, "total edges")
  TEST.LISTS.EQ(nrow(A), 20, "total nodes")
  TEST.LISTS.EQ(sum(colSums(A)!= 2), 0, "only two nodes per edge")

  rs = rowSums(A)
  i = which.max(rs)
  ## i should be the O-A pair
  ## Test most-connected pair.
  pair.i = pair.code.to.pair(rke$pc[i])
  Bi = (pair.i$donor=="O" && pair.i$patient=="A")
  TEST.BOOL(Bi, "testing most connected node")
}
test.remove.pairs = function(args) {
   n = args$n
   rke = rrke(n)
   rmv= sample(as.integer(n/2), 1)
   pair.ids = sample(1:n, rmv, replace=F)
   rke2 = remove.pairs(rke, pair.ids)
   # 1. test size of remainder
   TEST.LISTS.EQ(get.size(rke2), n-rmv, "size of remainder")
   # 2. test   new-edges =  old-edges - removed.edges
   rmv.edges = get.incident.edges(rke, pair.ids)
   TEST.LISTS.EQ(length(rke.edges(rke2)),  length(rke.edges(rke)) - rmv.edges, "new edges 1")
   
   kept.pairs = setdiff(1:n, pair.ids)
   kept.edges = get.internal.edges(rke, kept.pairs)
   # 3. test  new-edges = kept-edges
   TEST.LISTS.EQ(length(kept.edges),  length(rke.edges(rke2)), "new edges 2")
   
   for(type in c("U", "O", "R", "S")) {
     t.all = length( filter.pairs.by.type(rke, type) )
     t.removed = length(which( sapply(pair.ids, 
                                       function(i) pair.type(pair.code.to.pair(rke$pc[i])))==type))
     t.left = length( filter.pairs.by.type(rke2, type))
     TEST.LISTS.EQ(t.all, t.left + t.removed, sprintf(" %s pairs", type) )
   }
   return(T)
}
test.get.incident.edges = function(args) {
  n = args$n
  rke = rrke(n)
  pair.ids = sample(1:n, 10, replace=F)
  A = get.model.A(rke)
  edges = get.incident.edges(rke, pair.ids)
  Ai = A[pair.ids,]
  nedges.theor = length(which(colSums(Ai)>0))
  # 1. Ad-hoc = function output
  TEST.LISTS.EQ(nedges.theor, length(edges))
}
test.get.external.edges = function(args) {
  m = 3
  n=30
  rke.list = rrke.many(m=m, n=n, uniform.pra=T)
  rke.all = pool.rke(rke.list)
  hid=sample(3, 1)
  h.pairs = get.hospital.pairs(rke.all, hid)
  not.h.edges = get.external.edges(rke.all, h.pairs)
  
  m = max.matching(rke.all, remove.edges = not.h.edges)
  ## Test whether splitting the rke.all through  get.external 
  ## or through the list give the same matching (should be because they are the same graphs.)
  TEST.LISTS.EQ(m$matching$utility, 
                max.matching(rke.list[[hid]])$matching$utility, 
                "utility of 1 hospital")

}
## Tests for matching.R
test.max.matching = function(args) {
  for(filename in  list.files("tests",full.names=T,pattern="rke-res")) {
    load(filename)
    TEST.LISTS.EQ(res$util,  max.matching(res$rke)$matching$utility)
    print(sprintf("Checked against %s [OK]", filename) )
    rm(res)
  }
  return(T)
}
test.get.subgraph = function(args) {
  n = args$n
  rke = rrke(n)
  for(t in c("R","S")) {
    rke.sub   = get.subgraph(rke, type=t)
    types =  get.pair.types(rke.sub, rke.pairs(rke.sub))
    TEST.SETS.EQ(types, c(t), str=sprintf("Checking %s subgraph", t))
    
    TEST.LISTS.EQ(get.size(rke.sub), length(filter.pairs.by.type(rke, type=t)))
  }
  return(T)
}


### Tests for mechanism.R
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
test.rCM = function(args) {
  m = args$m 
  n = args$n
  rke.list = rrke.many(m=m, n=n, uniform.pra=T) 
  rke.all = pool.rke(rke.list)
  # 1. All truthful
  all.t = paste(rep("t", m), collapse="")
  U = rCM(rke.list, rke.all, all.t  )
  m.all = max.matching(rke.all)
  
  TEST.LISTS.EQ(m.all$matching$utility, sum(U), "all T-- total utility")
 
  ## 2. All deviate
  all.c = paste(rep("c", m), collapse="")
  U = rCM(rke.list, rke.all, all.c )
  m.ind = sapply(1:m, function(hid) max.matching(rke.list[[hid]])$matching$utility)
  TEST.LISTS.GEQ(U, m.ind, "rCM at least as good as selfish matching")
  
}
test.g.share = function(args) {
  n = 10
  k = 20
  z = rep(k, n)
  x = k * n
  TEST.LISTS.EQ( g.share(z, x),  rep(20, n), "all same")
  
  z[1] = 21
  TEST.LISTS.EQ( g.share(z, x),  c(21, rep(20, n-1)), "one more")
  
  z = c(1,1,1,1,100)
  x = 5
  TEST.LISTS.EQ( g.share(z, x),  rep(1, 5), "someone asking a lot")
  
  z = c(1,1,1,1,100)
  x = 104
  TEST.LISTS.EQ( g.share(z, x),  z, "someone asking a lot + can cover")
  
  z = c(1,1,1,1,100)
  x = 0
  TEST.LISTS.EQ( g.share(z, x),  rep(0, 5), "zero supply")
  
  z = c(0,0,0,0,0)
  x = 1000
  TEST.LISTS.EQ( g.share(z, x),  rep(0, 5), "zero demand")
  
  z = c(0,0,0,0,1000)
  x = 1000
  TEST.LISTS.EQ( g.share(z, x),  c(0,0,0,0,1000), "one guy")
  
  ## Check whether the +1 is distributed fairly
  z = c(2,2,2,2,2)
  x = 6
  i = sample(5, 1)
  reps= replicate(1000, {  which.max(g.share(z,x))==i})

  TEST.MEAN.EQ(reps, mu0=1/5, str="20% prob.")
  
  
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

test.Bonus.QS = function(args) {
  m = args$m
  n = args$n 
  
  rke.list = rrke.many(m=m,n=n,uniform.pra=T)
  rke.all = pool.rke(rke.list)
  QS = Bonus.QS(rke.all)
  
  
  ud.pcs = pair.codes.per.type("U")
  
  
  # Test the Q-part (#X-Y /hospital)
  for(h in 1:m) {
    rke.h = rke.list[[h]]
    TEST.LISTS.EQ(sum(QS$Q[[h]][-ud.pcs]), 0, "non-UD pairs are 0")
    
    for(i in ud.pcs) {
      pair.i = pair.code.to.pair(i)
      ud.theor = length(filter.pairs.by.donor.patient(rke.h,dtype=pair.i$donor, ptype=pair.i$patient))
      
      TEST.LISTS.EQ(ud.theor, QS$Q[[h]][i], str=sprintf(" code=%d", i))
    }
  }
  
  rm(list=c("rke.list", "rke.all"))
  ## Ad-hoc test here.
  load(file="tests/rkes-m3-n8.Rdata")
  rke.list = rke.many$rke.list
  rke.all = rke.many$rke.all
  QS = Bonus.QS(rke.all)
  
  TEST.LISTS.EQ(QS$Q[[1]][5], 2, str="H1 A-O")
  TEST.LISTS.EQ(QS$Q[[2]][5], 2, str="H2 A-O")
  
  TEST.LISTS.EQ(QS$Q[[1]][14], 0, str= "H1 AB-A")
  TEST.LISTS.EQ(QS$Q[[1]][15], 0, str="H1 AB-B")
  TEST.LISTS.EQ(QS$Q[[3]][13], 2, str="H3 AB-O")
  TEST.LISTS.EQ(QS$Q[[3]][9], 1, str="H3 B-O")
  TEST.LISTS.EQ(QS$Q[[2]][14], 0, str="H2 AB-A")
  
  
  x = c(1,0,0,0,0)
  z = rep(0, 5)
  TEST.LISTS.EQ( sapply(ud.pcs, function(i) length(QS$S[[1]][[i]])), x, str="UD matches in H1")
  TEST.LISTS.EQ( sapply(ud.pcs, function(i) length(QS$S[[2]][[i]])), z, str="UD matches in H2")
  TEST.LISTS.EQ( sapply(ud.pcs, function(i) length(QS$S[[3]][[i]])), z, str="UD matches in H3")
    
  
  return(T)
}
