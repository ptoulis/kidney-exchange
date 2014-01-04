## mechanisms_test.R
CHECK.rke.equal <- function(rke1, rke2) {
  CHECK_EQ(rke.size(rke1), rke.size(rke2))
  CHECK_SETEQ(rke1$pairs$pair.id, rke2$pairs$pair.id)
  CHECK_EQ(sum(rke1$edges$can.donate), sum(rke2$edges$can.donate))
  CHECK_SETEQ(rke1$pairs$pc, rke2$pairs$pc)
}

# Runs all tests for mechanisms.R
run.all.mechanisms.tests <- function(ntrials) {
  ok <- function() print("[OK]")
  for(i in 1:ntrials) {
    print(sprintf("#####   (%d/%d)   TESTs mechanisms.R    #####", i, ntrials))
    test.kpd.create()
    ok()
    test.play.strategies()
    ok()
    test.run.mechanism()
    ok()
    test.compute.ir.constraints()
    ok()
    test.compute.Rsubgraph.constraints()
    ok()
    test.OAB.matched.in.OUU()
  }
}

test.kpd.create <- function() {
  # Tests whether the KPD markets are created correctly.
  # Test 1: For the strategy profile "ttttt" check the reported
  #         RKE of some hospital is also the real RKE
  # Test 2: Test a similar condition when the hospital deviates.
  #
  pool = rrke.pool(m=6, n=15, uniform.pra=T)
  
  kpd.t = kpd.create(pool, "tttttt", include.3way=F)
  kpd.c = kpd.create(pool, "cccccc", include.3way=T)
  
  hid = sample(1:6, size=1)
  CHECK.rke.equal(kpd.t$reported.pool$rke.list[[hid]], pool$rke.list[[hid]])
  CHECK_SETEQ(subset(kpd.t$reported.pool$rke.all$pairs, hospital==hid)$pair.id,
              pool$rke.list[[hid]]$pairs$pair.id)  
  m = max.matching(pool$rke.list[[hid]], include.3way=T)
  rke.remainder = rke.remove.pairs(pool$rke.list[[hid]], rm.pair.ids=m$match$pair.id)
  CHECK_EQ(rke.size(rke.remainder), rke.size(kpd.c$reported.pool$rke.list[[hid]]))
}

test.play.strategies <- function() {
  # Check if hide ={} empty set when hospital is truthful
  # Check when hospital is deviating that |hide| = |max.matching|
  pool = rrke.pool(m=3, n=20, uniform.pra=T)
  CHECK_EXCEPTION({ play.strategies(pool$rke.list, "tt", include.3way=T) }, msg="Not correct #strategies")
  out = play.strategies(pool$rke.list, "ttt", include.3way=F)
  
  hid = sample(1:3, size=1)
  rke = pool$rke.list[[hid]]
  CHECK_SETEQ(rke.pair.ids(rke), out[[hid]]$report)
  CHECK_EQ(length(out[[hid]]$hid), 0, msg="Hospital is truthful")
  
  
  random.strategies = sample(c("t", "c", "c"))
  include.3way = sample(c(T,F), 1)
  # print(random.strategies)
  out = play.strategies(pool$rke.list, paste(random.strategies, collapse=""), include.3way=include.3way)
  hid = which(random.strategies=="c")[1]
  rke = pool$rke.list[[hid]]
  m = max.matching(rke, include.3way=include.3way, regular.matching=T)
  
  str = out[[hid]]
  CHECK_EQ(m$utility, length(str$hide), msg="|hide| = |max matching|")
  CHECK_SETEQ(rke.pair.ids(rke), union(str$report, str$hide))
  CHECK_DISJOINT(str$report, str$hide)
  # |matching| = 0 in the reported subgraph
  remainder = rke.keep.pairs(rke, pair.ids=str$report)
  CHECK_EQ(max.matching(remainder, include.3way=include.3way)$utility, 0)
}

test.run.mechanism <- function() {
  m = 2 + 2 * rpois(1, lambda=2)
  n = rpois(1, lambda=10)
 
  pool = rrke.pool(m=m, n=n, uniform.pra=F)
  include.3way = sample(c(T,F), 1)
  kpd = kpd.create(pool, strategy.str=paste(rep("c", m), collapse=""), include.3way=include.3way)
  
  mech = sample(c("rCM", "xCM", "Bonus"), 1)
  out = Run.Mechanism(kpd, mech, include.3way=include.3way)
  print(sprintf("TEST Run.Mechanism. m=%d hospitals, n=%d pairs, mech=%s, 3way=%s", m, n, mech, include.3way))
  mech.matching = do.call(mech, args=list(pool=kpd$reported.pool,
                                                      include.3way=include.3way))
  # 1. Test utilities
  CHECK_TRUE(abs(mech.matching$utility - out$mech.matching$utility) <= 3)
  
  internals = sapply(1:length(pool$rke.list), function(h) max.matching(pool$rke.list[[h]],
                                                                       include.3way=include.3way,
                                                                       regular.matching=T)$utility)
  CHECK_GE(get.matching.hospital.utilities(out$total.matching, m),
           internals)
}

test.compute.ir.constraints <- function() {
  # tests the IR constraints.
  # This is a data-frame
  #   pc   hospital internal.matches,    desc  pair.type
  pool = rrke.pool(m=rpois(1, lambda=6), n=rpois(1, lambda=20), uniform.pra=T)
  ir = compute.ir.constraints(pool, pair.types=c("S", "R"))
  CHECK_MEMBER(names(ir), c("pc", "hospital", "internal.matches", "desc", "pair.type"))
  print(sprintf("TEST compute.ir.constraints() -- m=%d, n=%d, Checking %d constraints. Avg=%.2f constraints", 
                length(pool$rke.list), rke.size(pool$rke.list[[1]]), nrow(ir),
                mean(ir$internal.matches)))
  for(i in 1:nrow(ir)) {
    constraint = ir[i, ]
    arg.desc = constraint$desc
    type = constraint$pair.type
    h = constraint$hospital
    
    h.pair.ids = subset(pool$rke.all$pairs, hospital==h & pair.type==type)$pair.id
    rke.h = rke.keep.pairs(pool$rke.list[[h]], pair.ids=h.pair.ids)
    
    m = max.matching(rke.h)
    nmatched = nrow(subset(as.data.frame(m$match), desc==arg.desc))
    # print(subset(m$match, desc==arg.desc))
    # print(sprintf("Desc=%s  type=%s hospital=%d,  matched internally=%d", arg.desc, type, h, nmatched))
    CHECK_EQ(constraint$internal.matches, nmatched)
  }
  
}

test.compute.Rsubgraph.constraints <- function() {
  pool = rrke.pool(m=rpois(1, lambda=6), n=rpois(1, lambda=20), uniform.pra=T)
  ir = compute.ir.constraints(pool, pair.types=c("S", "R"))
  Rsub = compute.Rsubgraph.constraints(ir, pool)
  out = join(Rsub, ir, by=c("pc", "hospital"), type="inner")
  
  internal.index = which(names(out) == "internal.matches")
  CHECK_EQ(length(internal.index), 2)
  if(nrow(out) > 0) {
    CHECK_EQ(out[, internal.index[1]], out[, internal.index[2]])
    
    totalPc = out$unmatched + out$internal.matches
    out.totalPc = sapply(1:nrow(out), function(i) nrow(subset(pool$rke.all$pairs, hospital==out[i,]$hospital & pc==out[i,]$pc)))
    CHECK_EQ(out.totalPc, totalPc)
  }
  return(out)
  
  
}

test.OAB.matched.in.OUU <- function() {
  rke = rrke(80)
  x = OAB.matched.in.OUU(rke)
  print(sprintf("TEST OAB.matched.in.OUU() -- Total OUU matches =%d", length(c(x$OAB, x$U))/3))
  CHECK_EQ(length(x$U), 2 * length(x$OAB))
  CHECK_SETEQ(names(x), c("OAB", "U", "matching"))
  nOAB = nrow(subset(rke$pairs, desc=="O-AB"))
  CHECK_GE(nOAB, length(x$OAB))
  if(length(x$U) > 0) {
    CHECK_SETEQ(x$matching$match$pair.type, c("O", "U"))
    CHECK_SETEQ(x$matching$matched.cycles$type, c(3))
  }
}