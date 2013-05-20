## Matching functions. 
##########################################################################
##   Maximum matching ##
## Returns { gurobi =>     { objvalue,  x = [0,1,...] },
##          matching =>  {matched.ids, matched.edges...} 
## x_i = 1  only if  edge i  is included in the maximum matching.
## Shuffling the edges eliminates bias in edge matching.
library(gurobi)
##  Helper function to max.matching()
get.matched.ids <- function(model.A, edge.ids) {
  ids = c()
  for(edge in edge.ids) {
    edge.ids = which(model.A[, edge]==1)
    ids = c(ids, edge.ids)
  }
  return(unique(ids))
}

##   Maximum  2min / maximum matching.
##  Can return NA if time out.
max.matching <- function(rke, 
                         regular.matching=F,
                         IR.constraints=list(),
                         shuffle.edges=T,
                         remove.edges=c(),
                         timeLimit=120) {
  
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
    model.obj.coefficients[remove.edges] <- -0.5 
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
  
  ##  Seems to be much faster than the old params.
  params.new <- list(OutputFlag=0,
                     NodefileStart=0.4,
                     Cuts=3,
                     Presolve=1,
                     MIPFocus=1,
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
  result$matching$timeout = F
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