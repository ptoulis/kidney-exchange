# TODO(ptoulis): Documentation
library(gurobi)

map.gurobiResult <- function(gurobi.result, rke, cycles) {
  # Gets the Gurobi output and returns the subset of rke "pairs" object
  # of those that have been matched.
  #
  # Returns: A "pairs" object.
  if(gurobi.result$status=="TIME_LIMIT" | gurobi.result$status=="INF_OR_UNBD") {
    warning("Time limit or infinity.")
    empty.result = get.empty.result()
    empty.result$status = gurobi.result$status
    return(empty.result)
  }
  
  # The convention is:  x[2-cycles, 3-cycles]
  matched.cycle.ids = which(gurobi.result$x == 1)
  # print (matched.cycle.ids)
  matched.cycles = cycles[matched.cycle.ids, ]
  matched.ids <- c(matched.cycles$pair.id1,
                   matched.cycles$pair.id2, 
                   matched.cycles$pair.id3)
  matched.ids <- setdiff(matched.ids, c(0))  # remove 0's (id3=for 2-way xchange)
  not.matched.ids = setdiff(rke.pair.ids(rke), matched.ids)
  CHECK_TRUE(all(!duplicated(matched.ids)), msg="No duplicates in matched ids")
  CHECK_DISJOINT(matched.ids, not.matched.ids, msg="Either matched or not")
 
  # Checks
  CHECK_EQ(length(matched.ids), 2 * nrow(subset(matched.cycles, type == 2)) + 
                                3 * nrow(subset(matched.cycles, type == 3)))
 
  # Return final result.
  result = subset(rke$pairs, pair.id %in% matched.ids)
  return (result)
}

##   Maximum  2min / maximum matching.
##  Can return NA if time out.
max.matching <- function(rke, include.3way=F,
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
  Cycles = rke.cycles(rke, include.3way=include.3way)
  model.w <- Cycles$type
  model.rhs <- rep(1, num.pairs)
  model.sense <- rep("<=", num.pairs)
  # Build A matrix of model (not confused with adj matrix of rke)
  # A2 =  pairs x 2cycles   A2ij = 1 if pair i in j 2-way cycle
  # Similarly, A3ij = 1 if pair i is in j 3-way exchange.
  model.A <- rke.cycles.membership(rke, Cycles)
  CHECK_EQ(nrow(model.A), num.pairs)
  CHECK_EQ(ncol(model.A), nrow(Cycles))
  CHECK_MEMBER(unique(model.A), c(0,1))
    
  model <- list()
  model$A          <- model.A
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
  match.out = map.gurobiResult(gurobi.result, rke, Cycles)
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