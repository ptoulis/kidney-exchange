## mechanisms_test.R

test.kpd.create <- function() {
  pool = rrke.pool(m=6, n=15, uniform.pra=T)
  
  kpd.t = kpd.create(pool, "tttttt", include.3way=F)
  kpd.c = kpd.create(pool, "cccccc", include.3way=T)
  test.rke.equal <- function(rke1, rke2) {
    CHECK_EQ(rke.size(rke1), rke.size(rke2))
    CHECK_SETEQ(rke1$pairs$pair.id, rke2$pairs$pair.id)
    CHECK_EQ(sum(rke1$edges$can.donate), sum(rke2$edges$can.donate))
    CHECK_SETEQ(rke1$pairs$pc, rke2$pairs$pc)
  }
  
  hid = sample(1:6, size=1)
  test.rke.equal(kpd.t$reported.pool$rke.list[[hid]], pool$rke.list[[hid]])
  
  m = max.matching(pool$rke.list[[hid]], include.3way=T)
  rke.remainder = rke.remove.pairs(pool$rke.list[[hid]], rm.pair.ids=m$match$pair.id)
  CHECK_EQ(rke.size(rke.remainder), rke.size(kpd.c$reported.pool$rke.list[[hid]]))
}

test.compute.ir.constraints <- function() {
  # tests the IR constraints.
  # This is a data-frame
  #   pc   hospital internal.matches,    desc  pair.type
  pool = rrke.pool(m=rpois(1, lambda=6), n=rpois(1, lambda=20), uniform.pra=T)
  ir = compute.ir.constraints(pool, pair.types=c("S", "R"))
  CHECK_MEMBER(names(ir), c("pc", "hospital", "internal.matches", "desc", "pair.type"))
  print(sprintf("m=%d, n=%d, Checking %d constraints", 
                length(pool$rke.list), rke.size(pool$rke.list[[1]]), nrow(ir)))
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
  CHECK_EQ(out[, internal.index[1]], out[, internal.index[2]])
  
  totalPc = out$unmatched + out$internal.matches
  out.totalPc = sapply(1:nrow(out), function(i) nrow(subset(pool$rke.all$pairs, hospital==out[i,]$hospital & pc==out[i,]$pc)))
  CHECK_EQ(out.totalPc, totalPc)
  
  return(out)
  
  
}