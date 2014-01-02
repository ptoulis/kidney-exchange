## mechanisms_test.R
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