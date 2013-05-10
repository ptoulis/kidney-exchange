# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges

## Jan 2013, new version of rke.R
rm(list=ls())
library(gurobi)
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
              compact=matrix(0, nrow=0, ncol=0))
   return(obj)
}
##  We allow a compact representation of the B, P matrices
to.compact <- function(P,B) {
    n = nrow(P)
    comp = matrix(NA, nrow=n, ncol=n)
    comp =    lower.tri(P) * P +upper.tri(B,diag=T) * B
    return(comp)
}
##  Get the P matrix only
get.P <- function(rke) {
    ## Take the lower part
    P = lower.tri(rke$compact) * rke$compact
    return(P+t(P))
}
##  Get the B matrix only.
get.B <- function(rke) {
    B = upper.tri(rke$compact) * rke$compact
    return(B+t(B)+ diag(diag(rke$compact)))
}







## Samples a rrke object.
rrke <- function(n, sample.pra.fn = rpra, verbose=F) {
    
    # 1. Sample the pairs
    pairs.obj = rpairs.new(n, sample.pra.fn, verbose)
    
    pair.codes = pairs.obj$codes
    pras = pairs.obj$pras
    
    bin.PRA.matrix =sample.bin.pra.matrix(pras, pras,same.hospital=T, verbose)
    bin.B.matrix =  get.bin.blood.matrix(pair.codes, verbose)
    
    # 2. compact form
    Compact.Matrix =to.compact(bin.PRA.matrix, bin.B.matrix)
    
    ## Define object to return.
    obj = list()
    obj$pc = pair.codes
    obj$compact = Compact.Matrix
    obj$pras = pras
    
    
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
rrke.many <- function(k=3, n=60) {
    x = list()
    for(i in sample(1:k)) 
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
    }
    #print(ranges)
    P.all = sample.bin.pra.matrix(rke.all$pras, rke.all$pras,same.hospital=T)
    B.all = get.bin.blood.matrix(rke.all$pc)
    for(i in 1:k) {
        s = ranges[i,1]
        t = ranges[i,2]
        P.all[s:t, s:t] = get.P(rke.list[[i]])
    }
    
    rke.all$compact=to.compact(P.all, B.all)
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
remove.pairs <- function(rke, pairs) {
   if(length(pairs)==0)
       return(rke)
   
   rke.new = list()
   rke.new$pras = rke$pras[-pairs]
   rke.new$pc = rke$pc[-pairs]
   rke.new$compact = rke$compact[-pairs, -pairs]
   
   return(rke.new)
}

## Returns the valus of the attribute "attr" for all pairs.
##  e.g. get.pairs.attribute(rke, ")
get.pairs.attribute <- function(rke, attr) {
    if(attr !="type")
        stop("Don't know how to do this. rke2.R line 144")
    
    return(sapply(1:get.size(rke), function(i) pair.type(pair.code.to.pair(rke$pc[i]))))
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
get.incident.nodes = function(A, edges) {
    ret = c()
    for(e in edges) {
        ret = c(ret, which(A[,e]==1))
    }
    return(unique(ret))
}
##   new function
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
        pairMatch = (type1 %in% t1 && type2 %in% t2) ||
                     (type1 %in%  t2 &&type2 %in% t1)
        return(pairMatch)
    })
   return(which(membership==T))
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
get.model.A <- function(rke, CAP=Inf) {
    ##  the adjacency matrix. 
    Adj = get.P(rke) * get.B(rke)
    Adj = upper.tri(Adj) * Adj ## get the upper triangular
    edges = c(0, cumsum(rowSums(Adj)) )  ##  n+1   elements.
    n = get.size(rke)
    ## Now make this a model matrix 
    K = sum(Adj)  ## total no. of edges
      
    A = matrix(0, nrow=n, ncol=K)
    kill.edges=  c()
    
    for(i in 1:n) {
        neighbors=  which(Adj[i, ]==1)
        right = edges[i+1]
        left = edges[i]+1

        if(right>=left) {
            #edges.to.add = right-left+1
            #edges.to.add = min(CAP, edges.to.add)
            #neighbors = sample(neighbors, edges.to.add)
            edges.to.add = left:right
            if(length(edges.to.add)>CAP) {
                n.kill = length(edges.to.add)- CAP
                kill.edges=  c(kill.edges, sample(edges.to.add, n.kill))
            }
            for(j in 1:length(edges.to.add))  {
                index= c(i, neighbors[j])
                A[index, edges.to.add[j] ]=1
            }
        }
    }
    if(CAP<Inf) {
        print(sprintf("Edges Before %d. After %d. Reduction  %.2f%%", 
                      ncol(A),
                      ncol(A)-length(kill.edges), 100* (length(kill.edges)/ncol(A))))
        warning("Capping is experimental")
    }
    if(length(kill.edges)>0)
        A = A[,-kill.edges]
    return(A)
    
}


##########################################################################
##   Maximum matching ##
## Returns { gurobi =>     { objvalue,  x = [0,1,...] },
##          matching =>  {matched.ids, matched.edges...} 
## x_i = 1  only if  edge i  is included in the maximum matching.
## Shuffling the edges eliminates bias in edge matching.

##   Maximum  2min / maximum matching.
##  Can return NA if time out.
max.matching <- function(rke, 
                        regular.matching=F,
                        IR.constraints=matrix(0,nrow=0, ncol=0),
                        shuffle.edges=T,
                        remove.edges=c(),
                        CAP=Inf, use.cutoff=F, 
                        timeLimit=120) {
    
    ## Size of RKE
    n = get.size(rke)
    ###   1.   Get the model matrix. 
    model.A = get.model.A(rke, CAP=CAP)
    ## Total no. of edges
     
    K = ncol(model.A) 
   
    get.empty.result <- function() {
        return(list(gurobi=list(objval=0, x=c()),
                    matching=list(matched.edges=c(),
                                  matched.ids=c(),
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
    
    
    ## Required regular matching. Put more weights on O-U edges
    if(regular.matching) {
        OUedges = filter.edges.by.type(rke, "O", "U")
        if(length(OUedges)>0)
            model.obj.coefficients[OUedges] = 2
    }
    if(length(remove.edges)>0) {
        model.obj.coefficients[remove.edges]=0
    }
    ###   If IR constraints  (used by xCM matching)
    ## IR.constraints = k x 2  matrix
    ## e.x. X[3, 1] = 5 means  for hospital 3 match at least as 5 "R" pairs.
    if(nrow(IR.constraints)>0) {
        ### Impose IR constraints only on R and S's
        if(length(rke$hospital)==0)
            stop("Cannot run constraints with < 2 hospitals")
        if(nrow(IR.constraints)!= length(unique(rke$hospital)))
            stop("Error. Row of IR constraints should = RKE hospitals")
        hospitals = 1:nrow(IR.constraints)
        h.size = length(hospitals)
        ## Returns 0,1,2  = # of pairs of <type> in <edge.id> for hospital <hid>
        get.edge.coefficient = function(edge.id, hospital.id, type) {
            edge.nodes = which(model.A[,edge.id]==1)
            s= 0
            for(node.j in 1:2) {
                ## get pair of edge
                pair = edge.nodes[node.j]
                ## get hospital of pair
                h= rke$hospital[pair]
                ## if match hospital AND match type add +1
                if(h==hospital.id && pair.type(pair.code.to.pair(rke$pc[pair]))==type)
                    s = s+1
            }
            return(s)
        }
        ## Constraints matrix  Hospitals x  K   
        ## For every hospital match at least as many internal matches of R,S
        ## A2.Rij = {0,1,2} whether edge j matches 0,1,2 R pairs for hospital i 
        A2.R = matrix(0, nrow=h.size, ncol=K)
        A2.S =  matrix(0, nrow=h.size, ncol=K)
        
        for(hid in hospitals) {
            A2.R[hid,]=sapply(1:K, function(edge.i) 
                get.edge.coefficient(edge.i, hid, "R" ))
            A2.S[hid,]=sapply(1:K, function(edge.i) 
                get.edge.coefficient(edge.i, hid, "S" ))
        }
        
        ## 2. Expand the constraints:
        ##    Match (R or S) at least as much as defined in IR.constraints
        model.A = rbind(model.A, A2.R)
        model.rhs = c(model.rhs, IR.constraints[,1])
        model.sense = c(model.sense, rep(">=", h.size))
        model.A = rbind(model.A, A2.S)
        model.rhs = c(model.rhs, IR.constraints[,2])
        model.sense = c(model.sense, rep(">=", h.size))
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
                       Threads=1,
                       Cuts=3,
                       Presolve=1,
                       TimeLimit=timeLimit)
    if(use.cutoff && CAP==Inf) {
        #m2= max.matching(rke, CAP=1)
        #cutoff = m2$matching$utility
        
        theor.cutoff= 0.5 * n - sqrt(n) - 2
        cutoff=theor.cutoff/3
        print(sprintf("Calculated cutoff= %.1f. Theoretical=%.1f", cutoff, theor.cutoff))
        params.new$Cutoff= cutoff
        
    }
    
      
    gurobi.result <- gurobi(model, params.new)
    if(gurobi.result$status=="TIME_LIMIT")
    {
        empty.result = list()
        empty.result$matching = list(matched.edges=c(), 
                                     matched.ids=c(), 
                                     not.matched.ids = c(),
                                     utility=0)
        empty.result$gurobi = list()
        return(empty.result)
    }
    
    old.x = rep(0,K)
    ## Fix shuffling
    for(j in 1:K) {
        old.edge = edges.shuffled[j]
        old.x[old.edge] = gurobi.result$x[j]
    }
    gurobi.result$x = old.x
    ##########
    matched.edges = which(gurobi.result$x==1)
    matched.ids = get.matched.ids(model.A, matched.edges)
    original.ids = 1:get.size(rke)
    not.matched.ids =  sort(setdiff(original.ids, matched.ids)) 
    
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