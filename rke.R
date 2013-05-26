# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges
## Jan 2013, new version of rke.R
rm(list=ls())
source("lib.R")
source("testing.R")
source("tables.R")
##      The main data structure.
##        rke2 = KE pool object (one hospital)
##               {  pc=> [1,2,1,1,4,13,...]    # pair codes
##                  pras=>[0.2, 0.3, 0.1, ...]
##                  B = matrix
##                  P = matrix
################################################################
#  Return an empty KE object
##   Defined by pc = vector of pair codes
##              pras = vector of PRA sensitivity
##              P = binary matrix of crossmatch
##              B = binary matrix of blood-type compatibility
empty.rke <- function() {
   obj = list(pc=c(),  
              pras=c(),
              P=matrix(0, nrow=0, ncol=0),
              B=matrix(0, nrow=0, ncol=0),
              hospital=c())
   return(obj)
}

## Samples a rrke object.
rrke <- function(n, 
                 uniform.pra = T,
                 blood.type.distr = list(O=0.5, A=0.3, B=0.15, AB=0.05),
                 verbose=F) {
    if(n==0) return(empty.rke() )
    # 1. Sample the pairs
    ##   each one has a code and a PRA
    pairs.obj = rpairs(n, 
                       uniform.pra, 
                       blood.type.distr,
                       verbose)
    
    pair.codes = pairs.obj$codes
    pras = pairs.obj$pras
    
    bin.PRA.matrix =rpra.matrix(pras, pras,same.hospital=T, verbose)
    bin.B.matrix =  get.bin.blood.matrix(pair.codes, verbose)
    
     ## Define object to return.
    obj = list()
    obj$pc = pair.codes
    obj$P = bin.PRA.matrix
    obj$B = bin.B.matrix
    obj$pras = pras
    obj$uniform.pra = uniform.pra
    ###  Checks 
    #print(sum(get.P(obj)==bin.PRA.matrix) == n^2)
    #print(sum(get.B(obj)==bin.B.matrix) == n^2)
    
    ## Clean-up
    rm(bin.B.matrix)
    rm(bin.PRA.matrix)
    rm(pras)
    rm(pair.codes)
    return(obj)
}




#  k = # hospitals
#  n = # pairs/hospital
rrke.many <- function(m=3, n=60, uniform.pra) {
  
  x = list()
  for(i in sample(1:m)) 
    x[[i]] = rrke(n,uniform.pra=uniform.pra)
  return(x)
}

## Total number of pairs in the graph
get.size <- function(rke) {
  return( length(rke$pc) )
}

##   Merge different RKE objects into one.
##   Assume all hospitals are non-empty
pool.rke <- function(rke.list) {
  
  rke.all = empty.rke()
  k  = length(rke.list)
  ranges = matrix(0, nrow=k, ncol=2)
  ranges[,2] = cumsum(sapply(1:k, function(i) get.size(rke.list[[i]])))
  ranges[,1] = c(1, ranges[-k,2]+1)
  for(i in 1:k) {
    rke.all$pc = c(rke.all$pc, rke.list[[i]]$pc)
    rke.all$pras =c(rke.all$pras, rke.list[[i]]$pras)
    rke.all$uniform.pra[i] = rke.list[[i]]$uniform.pra
  }
  #print(ranges)
  P.all = rpra.matrix(rke.all$pras, rke.all$pras,same.hospital=T)
  B.all = get.bin.blood.matrix(rke.all$pc)
  for(i in 1:k) {
    s = ranges[i,1]
    t = ranges[i,2]
    P.all[s:t, s:t] = rke.list[[i]]$P
  }
  
  rke.all$P  = P.all
  rke.all$B = B.all
  rke.all$hospital = c()
  for(j in 1:k) {
    rke.all$hospital = c(rke.all$hospital, rep(j, get.size(rke.list[[j]])))
  }
  # clean up
  rm(list=c("P.all", "B.all"))
  rke.all$map.back = function(new.ids) {
    ## should return the old hospital id
    ret = c()
    for(id in new.ids) {
      if(id > get.size(rke.all) || id <0)
        stop("invalid new.id in map.back()")
      hosp.id = rke.all$hospital[id]
      old.id = id - min(which(rke.all$hospital==hosp.id))+1
      ret = c(ret, old.id)
    }
    return( ret )
  }
  rke.all$map.fwd = function(old.id, hid) {
    return( min(which(rke.all$hospital==hid) ) + old.id-1 )
  }
  return(rke.all)
}

##  Given an RKE  (1) subtract (2) return remainder.
# Used to represent deviation strategies ("hide")
remove.pairs <- function(rke, pair.ids) {
  
  all.pairs = rke.pairs(rke)
  if(! is.subset(all.pairs, pair.ids))
    stop("Pair ids to remove should be a subset!")
  
  ## First create the map-back function.
  count = 0
  map.ids = c()
  for(i in 1:get.size(rke)) {
    if(i %in% pair.ids) {
      map.ids[i] = 0
    } else {
      count = count+1
      map.ids[i] = count
    }
  }
  ## map.ids= [1,2,0,0,0,3,0,0,4,5,0,0,6]
  
  map.back = function(new.ids) {
    ## should return the old hospital id
    ret = c()
    for(id in new.ids) {
      if(id > get.size(rke.new) || id <0)
        stop("invalid new.id in map.back()")
      ret = c(ret, which(map.ids==id))
    }
    return( ret )
  }
  
  ## define rke.new
  rke.new = list()
  
  if(equal.sets(pair.ids, all.pairs))
  { 
    rke.new = empty.rke() 
  } else if(length(pair.ids)==0) {
    rke.new = rke
  } else { 
    rke.new$pras = rke$pras[-pair.ids]
    rke.new$pc = rke$pc[-pair.ids]
    rke.new$P = rke$P[-pair.ids, -pair.ids]
    rke.new$B = rke$B[-pair.ids, -pair.ids]
    ## this is necessary. R messes up splicing! 
    ## if you take out n-1 of the dimensions it will return a vector instead of a matrix!!! argh!
    if(length(pair.ids) == length(all.pairs)-1) {
      rke.new$P = matrix(rke.new$P, nrow=1, ncol=1)
      rke.new$B = matrix(rke.new$B, nrow=1, ncol=1)
    }
    rke.new$uniform.pra = rke$uniform.pra
    if("hospital" %in% names(rke))
      rke.new$hospital = rke$hospital[-pair.ids]
  }
  
  rke.new$map.back = map.back  
  return(rke.new)
}
keep.pairs = function(rke, pair.ids) {
  warning("keep.pairs()  not unit-tested")
  all.pairs = rke.pairs(rke)
  rm.pairs = setdiff(all.pairs, pair.ids)
  return(remove.pairs(rke, rm.pairs))
}

## Returns the valus of the attribute "attr" for all pairs.
##  e.g. get.pairs.attribute(rke, ")
get.pairs.attribute <- function(rke, attr) {
  if(attr !="type")
    stop("Don't know how to do this. rke2.R line 144")
  
  return(sapply(1:get.size(rke), function(i) pair.type(pair.code.to.pair(rke$pc[i]))))
}

##  Returns a sub-RKE object of only type=OD/UD, R, S 
get.subgraph <- function(rke, type) {

  if(type=="S") {
    pairs.rmv = union(filter.pairs.by.type(rke, "R"), 
                      union(filter.pairs.by.type(rke, "O"), filter.pairs.by.type(rke,"U")))
    rke.new = remove.pairs(rke, pairs.rmv)
    return(rke.new)            
    
  }
  if(type=="R") {
    pairs.rmv = union(filter.pairs.by.type(rke, "S"), 
                      union(filter.pairs.by.type(rke, "O"), filter.pairs.by.type(rke,"U")))
    rke.new = remove.pairs(rke, pairs.rmv)
    return(rke.new)   
  }
  if(type=="O/U") {
    pairs.rmv = union(filter.pairs.by.type(rke, "R"), filter.pairs.by.type(rke,"S"))
    rke.new = remove.pairs(rke, pairs.rmv)
    return(rke.new)   
  }
  stop("Wrong type requested in get.subgraph")
}
get.incident.nodes = function(rke, edges) {
  ret = c()
  A = get.model.A(rke)
  if(max(edges)> length(rke.edges(rke)))
    stop("Invalid edges. ")
  for(e in edges) {
    ret = c(ret, which(A[,e]==1))
  }
  return(unique(ret))
}
get.incident.edges = function(rke, pair.ids) {
  
  A = get.model.A(rke)
  if(length(pair.ids) == 0 || ncol(A)==0) return(c())
  Ai = matrix( A[pair.ids,], ncol=ncol(A) )
  if(length(pair.ids)==1)
    return(which(Ai>0))
  
  return(which(apply(Ai, 2, sum)>0))
}
get.nonincident.edges = function(rke, pair.ids) {
  return(setdiff(rke.edges(rke), get.incident.edges(rke, pair.ids)))
}
get.internal.edges = function(rke, pair.ids) {
  pair.ids = c(pair.ids)
  if(length(pair.ids)<2) return(c())
  A = get.model.A(rke)
  if(ncol(A)==0) return(c())
  Ai = matrix( A[pair.ids,], ncol=ncol(A) )
  
  return(which(apply(Ai, 2, sum)==2))
}
get.external.edges = function(rke, pair.ids) {
  z = get.internal.edges(rke, pair.ids)
  if(length(z)==0)
    return(rke.edges(rke) )
  
  return(setdiff(rke.edges(rke), z))
}
##   new function
rke.edges = function(rke) {
  A = get.model.A(rke)
  if(ncol(A)==0) return(c())
  return(1:ncol(A))
}
rke.pairs = function(rke) {
  n= get.size(rke)
  if(n==0) return(c())
  return(1:n)
}

get.pair.types = function(rke, pair.ids) {
  warnings("get.pair.types not unit-tested")
  pair.ids = c(pair.ids)  ## avoid some stupid bugs
  if(length(pair.ids)==0) return(c() )
  return( get.pairs.attribute(rke, attr="type")[pair.ids])
}

##  FILTERS 
##  for pairs and edges.
filter.edges.by.type <- function(rke, t1, t2) {
  #print(sprintf("t1=%s t2=%s", t1, t2))
  A = get.model.A(rke)
  
  if(t1=="*") 
    t1 = c("O", "U", "S", "R")
  else t1 = c(t1)
  
  if(t2=="*") 
    t2 = c("O", "U", "S", "R")
  else t2=c(t2)
  
  K = ncol(A)
  if(K==0) return(c())
  
  membership = sapply(1:K, function(e) {
    ids = which(A[,e]==1)
    pair1=  pair.code.to.pair(rke$pc[ids[1]])
    pair2=  pair.code.to.pair(rke$pc[ids[2]])
    
    type1 = pair.type(pair1)
    type2 = pair.type(pair2)
    #print(sprintf("Edge %d t1, t2=%s, %s", e, type1, type2))
    #print(type2)
    match = (type1 %in% t1 && type2 %in% t2) ||
      (type1 %in%  t2 &&type2 %in% t1)
    return(match)
  })
  return(which(membership==T))
}
filter.edges.by.donor.patient <- function(rke, dt, pt) {
  #print(sprintf("t1=%s t2=%s", t1, t2))
  A = get.model.A(rke)
  
  if(dt=="*") 
    dt = c("O", "A", "B", "AB")
  else dt = c(dt)
  
  if(pt=="*") 
    pt = c("O", "A", "B", "AB")
  else pdt = c(pt)
  
  K = ncol(A)
  if(K==0) return(c())
  
  membership = sapply(1:K, function(e) {
    ## find incident nodes to edge
    edge.ids = which(A[,e]==1)
    pair1=  pair.code.to.pair(rke$pc[edge.ids[1]])
    pair2=  pair.code.to.pair(rke$pc[edge.ids[2]])
    
    pairMatch = (pair1$donor %in% dt && pair1$patient %in% pt) ||
      (pair2$donor %in% dt && pair2$patient %in% pt)
    return(pairMatch)
  })
  return(which(membership==T))
}
filter.out.edges.by.type <- function(rke, t1, t2) {
  edges = filter.edges.by.type(rke, t1, t2)
  return(setdiff(rke.edges(rke), edges))
}
##  Returns only those pairs of specific donor-patient types.
filter.pairs.by.donor.patient <- function(rke, dtype, ptype) {
  if(get.size(rke)==0) return(c() )
  if(dtype=="*") 
    dtype = c("O", "A", "B", "AB")
  else dtype = c(dtype)
  
  if(ptype=="*") 
    ptype = c("O", "A", "B", "AB")
  else ptype=c(ptype)
  
  n = get.size(rke)
  if(n==0) return(c())
  membership = sapply(1:n, function(i) {
    pair =  pair.code.to.pair(rke$pc[i])
    donor.type = pair$donor
    pat.type = pair$patient
    pairMatch = (donor.type %in% dtype && pat.type %in% ptype)
    return(pairMatch)
  })
  return(which(membership==T))
}
filter.pairs.by.type <- function(rke, type) {
  if(get.size(rke)==0) return(c()) 
  return(which(get.pairs.attribute(rke, "type")==type))
}

get.hospital.pairs <- function(rke.all, h.ids) {
  if(! "hospital" %in% names(rke.all))
    stop("Input is not a pooled RKE object")
  if(length(h.ids)==1)
    return(which(rke.all$hospital == h.ids))
  return(which(rke.all$hospital %in% h.ids))
}
get.model.A <- function(rke) {
  ##  the adjacency matrix. 
  Adj = rke$P * rke$B
  Adj = upper.tri(Adj) * Adj
  
  all.pairs = rke.pairs(rke)
  
  ## The model matrix
  ## TO-DO(ptoulis):  Probably not very efficient?
  ## total number of edges
  K = sum(Adj)
  A = matrix(0, nrow=length(all.pairs), ncol=K)
  edge.counter = 0
  for(i in all.pairs) {
    neighbors =  sort(which(Adj[i, ]==1))
    for(j in neighbors) {
      edge.counter = edge.counter+1
      ##   Add the pair, only if it has not been entered before.
      A[c(i,j), edge.counter] <- 1
    }
  }
  TEST.LISTS.EQ(edge.counter, K)
  return(A)
}

