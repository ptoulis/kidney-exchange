# TODO(ptoulis): Better documentation here.
# TODO(ptoulis): Implementation for directed graph.
library(gurobi)

empty.match.result <- function() {
  return(list(matched.edges=c(),
              matched.pairs=c(),
              utility=0,
              not.matched.pairs=c()))
}

map.gurobiResult <- function(gurobi.result) {
  if(gurobi.result$status=="TIME_LIMIT" | gurobi.result$status=="INF_OR_UNBD") {
    warning("Time limit or infinity.")
    empty.result = get.empty.result()
    empty.result$status = gurobi.result$status
    return(empty.result)
  }
  
  which(gurobi.result$x==1)
  return(gurobi.result$x)
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
  return(result)
}

##   Maximum  2min / maximum matching.
##  Can return NA if time out.
max.matching <- function(rke, enable.3way=F,
                         IR.constraints=list(),
                         timeLimit=120) {
  CHECK_rke(rke)
  num.pairs = rke.size(rke)
  num.edges = length(rke.edge.ids(rke))
  if (num.edges == 0) {
    warning("Empty RKE object")
    return (empty.match.result())
  }
  # Define Gurobi model 
  # Gurobi defines the problem as: 
  # A * x   <sense>   rhs   ,  sense in {"<=", ">="}
  #
  # Follow cycle-formulation.
  # matrix of cycles.
  Cycles2way = rke.2way.cycles(rke)
  Cycles3way = matrix(0, nrow=0, ncol=3)
  if (enable.3way) {
    Cycles3way = rke.3way.cycles(rke)
  }
  
  num.2wayCycles = nrow(Cycles2way)
  num.3wayCycles = nrow(Cycles3way)
  model.size = num.2wayCycles + num.3wayCycles
  model.w <- c(rep(2, num.2wayCycles), rep(3, num.3wayCycles))
  model.rhs        <- rep(1, num.pairs)
  model.sense      <- rep("<=", num.pairs)
  # Build A matrix of model (not confused with adj matrix of rke)
  # A2 =  pairs x 2cycles   A2ij = 1 if pair i in j 2-way cycle
  # Similarly, A3ij = 1 if pair i is in j 3-way exchange.
  model.A2cycle <- matrix(NA, nrow=num.pairs, ncol=0)
  model.A3cycle <- matrix(NA, nrow=num.pairs, ncol=0)
  pair.ids = rke.pair.ids(rke)
  model.A2cycle = t(sapply(pair.ids, function(i) rowSums(Cycles2way == i)))
  if (enable.3way) 
    model.A3cycle <- t(sapply(pair.ids, function(i) rowSums(Cycles3way == i)))
  CHECK_MEMBER(unique(model.A2cycle), c(0,1))
  CHECK_MEMBER(unique(model.A3cycle), c(0,1))
  CHECK_EQ(ncol(model.A2cycle), num.2wayCycles, "Correct #col=#2-way cycles")
    
  model <- list()
  model$A          <- cbind(model.A2cycle, model.A3cycle)
  model$obj        <- model.w
  model$modelsense <- "max"
  model$rhs        <- model.rhs
  model$sense      <- model.sense
  model$vtype      <- rep('B', length(model.w))
  
  ##  Seems to be much faster than the old params.
  params.new <- list(OutputFlag=0,
                     NodefileStart=0.4,
                     Cuts=3,
                     Presolve=1,
                     MIPFocus=2,
                     TimeLimit=timeLimit)
  
  gurobi.result <- gurobi(model, params.new)
  match.out = map.gurobiResult(gurobi.result)
  if(all(is.element(gurobi.result$x, c(0,1)))) {
    return(match.out)
  } else {
      # Sometimes GUROBI gives results that are not in (0,1)
      # although we have set the variables to be binary (maybe i am missing something)
      # For that reason we call the function again. Empirically, there is no danger of 
     # infinite recursion although this is something that needs to be fixed.
      logwarn("Gurobi unstable output. Saving problematic RKE object AND retrying..")
      save(rke, file="debug/unstable.Rdata")
      return (max.matching(rke=rke, enable.3way=enable.3way,
                           IR.constraints=IR.constraints, 
                           timeLimit=timeLimit))
  }
}