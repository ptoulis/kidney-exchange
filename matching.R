# TODO(ptoulis): Better documentation here.
# TODO(ptoulis): Implementation for directed graph.
library(gurobi)

empty.match.result <- function() {
  return(list(matched.edges=c(),
              utility=0,
              not.matched.pairs=c()))
}

map.gurobiResult <- function(gurobi.result, rke, cycles.2way, cycles.3way) {
  if(gurobi.result$status=="TIME_LIMIT" | gurobi.result$status=="INF_OR_UNBD") {
    warning("Time limit or infinity.")
    empty.result = get.empty.result()
    empty.result$status = gurobi.result$status
    return(empty.result)
  }
  
  # The convention is:  x[2-cycles, 3-cycles]
  num.2cycles = nrow(cycles.2way)
  matched.cycles = which(gurobi.result$x == 1)
  cycles2.matched = matched.cycles[matched.cycles <= num.2cycles]
  cycles3.matched = setdiff(matched.cycles, cycles2.matched) - num.2cycles
  
  pairs.matched2 = as.vector(cycles.2way[cycles2.matched, ])
  pairs.matched3 = as.vector(cycles.3way[cycles3.matched, ])
  original.ids = rke.pair.ids(rke)
  matched.ids = c(pairs.matched2, pairs.matched3)
  not.matched.ids =  setdiff(original.ids, matched.ids)
  # Checks
  CHECK_DISJOINT(pairs.matched2, pairs.matched3, "pairs in 2-cycles != pairs in 3-cycles")
  CHECK_UNIQUE(pairs.matched2, msg="Unique 2-cycle matches")
  CHECK_UNIQUE(pairs.matched3, msg="Unique 3-cycle matches")
  CHECK_MEMBER(pairs.matched2, original.ids, msg="Valid matched ids")
  CHECK_MEMBER(pairs.matched3, original.ids, msg="Valid matched ids")
  CHECK_EQ(length(matched.ids), 2 * length(cycles2.matched) + 3 * length(cycles3.matched))
  # Return final result.
  result = empty.match.result()
  matched.edges <- c()
  for (c2 in cycles2.matched) {
    ids = cycles.2way[c2,]
    edge.id1 = rke.edge.id(rke, ids[1], ids[2])
    edge.id2 = rke.edge.id(rke, ids[2], ids[1])
    matched.edges <- c(matched.edges, c(edge.id1, edge.id2))
  }
  for (c3 in cycles3.matched) {
    ids = cycles.3way[c3,]
    edge.id1 = rke.edge.id(rke, ids[1], ids[2])
    edge.id2 = rke.edge.id(rke, ids[2], ids[3])
    edge.id3 = rke.edge.id(rke, ids[3], ids[1])
    matched.edges <- c(matched.edges, c(edge.id1, edge.id2, edge.id3))
  }
  result$matched.edges <- matched.edges
  result$matched.ids = matched.ids
  result$not.matched.ids = not.matched.ids
  result$utility = length(matched.ids)
  result$timeout = F
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
  if (nrow(Cycles2way) > 0)
    model.A2cycle = t(sapply(pair.ids, function(i) rowSums(Cycles2way == i)))
  if (enable.3way & nrow(Cycles3way) > 0) 
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
  match.out = map.gurobiResult(gurobi.result, rke, Cycles2way, Cycles3way)
  edge.index = which(rke$edges$edge.id %in% match.out$matched.edges)
  rke$edges$edge.color[edge.index] <- rep("red", length(match.out$matched.edges))
  plot.rke(rke)
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