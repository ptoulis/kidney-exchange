# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges
## Jan 2013, new version of rke.R
rm(list=ls())
source("lib.R")
##      The main data structure.
##        rke2 = KE pool object (one hospital)
##               {  pc=> [1,2,1,1,4,13,...]    # pair codes
##                  pras=>[0.2, 0.3, 0.1, ...]
##                  compact=> [nxn]  matrix
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
        x[[i]] = rrke(n)
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
    
    return(rke.all)
}

##  Given an RKE  (1) subtract (2) return remainder.
# Used to represent deviation strategies ("hide")
remove.pairs <- function(rke, pair.ids) {
  ### this is a delicate process, so we add some extra checks.
   if(length(pair.ids)==0)
       return(rke)
   all.pairs = rke.pairs(rke)
   if(! is.subset(all.pairs, pair.ids))
     stop("Pair ids to remove should be a subset!")
    
   if(equal.sets(pair.ids, all.pairs))
     return(empty.rke())

   rke.new = list()
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
   return(rke.new)
}
keep.pairs = function(rke, pair.ids) {
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
  if(length(pair.ids)==0) return(c())
  A = get.model.A(rke)
  Ai = A[pair.ids,]
  if(length(pair.ids)==1)
    return(which(sum(Ai)>0))
  
  return(which(apply(Ai, 2, sum)>0))
}
get.nonincident.edges = function(rke, pair.ids) {
  return(setdiff(rke.edges(rke), get.incident.edges(rke, pair.ids)))
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
  return(which(get.pairs.attribute(rke, "type")==type))
}

#
#
##    Maximum   Matching  code.
#
#
##  Helper function to max.matching()
get.matched.ids <- function(model.A, edge.ids) {
    ids = c()
    for(edge in edge.ids) {
        edge.ids = which(model.A[, edge]==1)
        ids = c(ids, edge.ids)
    }
    return(unique(ids))
}

##    Get the model matrix.    Need for max matching
##   TO-DO   Aij = 1  iff   edge j is incident on node i 
get.model.A <- function(rke) {
    ##  the adjacency matrix. 
    Adj = rke$P * rke$B
    all.pairs = rke.pairs(rke)
    n = length(all.pairs)
    if(sum(Adj)%%2 != 0) 
      stop("Adjacency matrix is not symmetric?")
  
    ## The model matrix
    ## TO-DO(ptoulis):  Probably not very efficient?
    A = matrix(0, nrow=n, ncol=0)
    for(i in all.pairs) {
        neighbors =  which(Adj[i, ]==1)
        for(j in neighbors[neighbors>i]) {          
          ##   Add the pair, only if it has not been entered before.
          A = cbind(A, rep(0, n))
          A[c(i,j), ncol(A)] <- 1
          
        }
    }
    return(A)
}
##########################################################################
##   Maximum matching ##
## Returns { gurobi =>     { objvalue,  x = [0,1,...] },
##          matching =>  {matched.ids, matched.edges...} 
## x_i = 1  only if  edge i  is included in the maximum matching.
## Shuffling the edges eliminates bias in edge matching.
library(gurobi)
##   Maximum  2min / maximum matching.
##  Can return NA if time out.
max.matching <- function(rke, 
                        regular.matching=F,
                        IR.constraints=list(),
                        shuffle.edges=T,
                        remove.edges=c(),
                        timeLimit=120) {
    warning("max.matching does not have a unit test")
    ## Size of RKE  (# pairs)
    n = get.size(rke)
    ###   1.   Get the model matrix. 
    model.A = get.model.A(rke)
    
    ## Total no. of edges
    K = ncol(model.A) 
   
    get.empty.result <- function() {
        return(list(gurobi=list(objval=0, x=c()),
                    matching=list(matched.edges=c(),
                                  matched.ids=c(),
                                  utility=0,
                                  not.matched.ids=c())))
    }
    if(K==0) {
        #print("No edges. Max matching is empty");
        return( get.empty.result() )
    }
    
    ## Need to create  A matrix :  Aij = 1 iff edge j is incident on node i
    
    ###   Define Gurobi model 
    ##  The Gurobi defines the problem as: 
    ##   A * x   <sense>   rhs   ,  sense in {"<=", ">="}
    ##  First n constraints = RKE nodes constraints
    ##  1.  No node should be double counted in matching
    model.obj.coefficients = rep(1,K)
    model.rhs        <- rep(1, n)
    model.sense      <- rep("<=",n)
    
    
    ## Required regular matching. Put more weights on O-U edges (almost-regular)
    if(regular.matching) {
        OUedges = filter.edges.by.type(rke, "O", "U")
        if(length(OUedges)>0)
            model.obj.coefficients[OUedges] = 2
    }
    ##  Remove the specified edges.
    if(length(remove.edges)>0) {
        model.obj.coefficients[remove.edges]=0
        if(length(remove.edges)==K)
          return(get.empty.result())
    }
    ###   If IR constraints  (used by xCM matching)
    ## IR.constraints = [][]  , i.e 
    ## e.x. X[[3]] = [1,0,0,0,2,.....]    1 x 16  vector
    ##  Vector has #pairs to match for each pair-code (1 to 16)
    if(length(IR.constraints)>0) {
        ### Impose IR constraints only on R and S's
        if(length(rke$hospital)==0)
            stop("Cannot run constraints with no hospitals")
      
        ## Returns 0,1,2  = # of pairs of <pair> in <edge.id> for hospital <hid>
        # pair = list(donor="A", patient="B")
        get.edge.coefficient = function(edge.id, hospital.id, pair.code) {
            edge.nodes = which(model.A[,edge.id]==1)
            s= 0
            for(pair.id in edge.nodes) {
                ## get hospital of pair
                h= rke$hospital[pair.id]
                ## if match hospital AND match type add +1
                if(h==hospital.id && rke$pc[pair.id] == pair.code)
                    s = s+1
            }
            return(s)
        }
        ## Build the constrained matrix = A'ij = 1  if  i-th contraint include edge j
        hospitals = 1:length(IR.constraints)  
        A.ConstrainedEdges = matrix(0, nrow=0, ncol=K)
        Rhs.ConstrainedEdges = matrix(0, nrow=0, ncol=1)
        ###   TO-DO(ptoulis): For loops +rbind makes things slow. Improve?
        for(hid in hospitals) {
          ## this is 1x16  vector of lower bounds for matches.
          ir.constraints.h = IR.constraints[[hid]]
          active.pcs = which(ir.constraints.h>0)
          
          for(pc in active.pcs) {
              pc.edge.coeffs = sapply(1:K, function(edge.i) get.edge.coefficient(edge.i, hid, pc))
              A.ConstrainedEdges = rbind(A.ConstrainedEdges, pc.edge.coeffs)
              Rhs.ConstrainedEdges = rbind(Rhs.ConstrainedEdges, ir.constraints.h[pc])
          }
         
        }
        ## 2. Expand the constraints:
        ##    Match at least as much as defined in IR.constraints
        model.A = rbind(model.A, A.ConstrainedEdges)
        model.rhs = c(model.rhs, Rhs.ConstrainedEdges)
        model.sense = c(model.sense, rep(">=", nrow(Rhs.ConstrainedEdges)))
    }### If   IR constraints
    
    ### Shuffle Edges (if required)
    edges.shuffled= 1:K
    if(shuffle.edges) 
        edges.shuffled = sample(sample(1:K)) 
    
    model <- list()
    model$A          <- as.matrix(model.A[,edges.shuffled])
    model$obj        <- model.obj.coefficients[edges.shuffled]
    model$modelsense <- "max"
    model$rhs        <- model.rhs
    model$sense      <- model.sense
    model$vtype      <- rep('B', K)
    
    params.def  <- list(OutputFlag=0)
    params.old<- list(OutputFlag=0,
                       NodefileStart=0.4,
                       Threads=1,
                       Cuts=3,
                       Presolve=1)
    
    ##  Seems to be much faster than the old params.
    cutoff = -Inf
    params.new <- list(OutputFlag=0,
                       NodefileStart=0.4,
                       Cuts=3,
                       Presolve=1,
                       MIPFocus=2,
                       TimeLimit=timeLimit)
  
    gurobi.result <- gurobi(model, params.new)
    if(gurobi.result$status=="TIME_LIMIT")
    {
        warning("TIME LIMIT hit in Gurobi run.")
        empty.result = get.empty.result()
        empty.result$TIMEOUT=T
        return(empty.result)
    }
    
    old.x = rep(0,K)
    ## Fix shuffling
    for(j in 1:K) {
        old.edge = edges.shuffled[j]
        old.x[old.edge] = gurobi.result$x[j]
    }
    gurobi.result$x = old.x
    ##########    Preparing the result
    matched.edges = which(gurobi.result$x==1)
    matched.ids = my.sort(get.matched.ids(get.model.A(rke), matched.edges) )
    original.ids = rke.pairs(rke)
    not.matched.ids =  my.sort(setdiff(original.ids, matched.ids)) 
    
    result = get.empty.result()
    result$matching$matched.edges= matched.edges
    result$matching$matched.ids = matched.ids###    TO -DO   count the matched. idsmatched.ids
    result$matching$not.matched.ids = not.matched.ids
    result$matching$utility = length(matched.ids)
    result$matching$timeout = T
    result$gurobi = gurobi.result
    #########
    
    if(  prod(unique(sort(result$gurobi$x)) %in% c(0,1))==1)
    {
        return(result)
    }   else {
        warning("Gurobi unstable output. Saving problematic RKE object AND retrying..")
        save(rke, file="debug/unstable.Rdata")
        #stop("Print submit this file to ptoulis@fas.harvard.edu")
        return (max.matching(rke, 
                             regular.matching,
                             IR.constraints,
                             shuffle.edges))
    }
    
}