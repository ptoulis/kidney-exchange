# TODO(ptoulis): Documentation
library(gurobi)

# A matching represents a (donor, patient) matching
#   match = a "pairs" object
#   utility = #transplants performed
#   status = "OK" if everything is fine
CHECK_matching <- function(matching) {
  # Checks whether a matching has the valid format.
  CHECK_MEMBER(c("match", "utility", "status", "matched.cycles"), names(matching), msg="Matching members")
  CHECK_pairs(matching$match)
  CHECK_GE(matching$utility, 0, "Utility is >=0")
  CHECK_TRUE(is.character(matching$status), msg="Status should be a string")
}

get.matching.from.ids <- function(matched.ids, rke) {
  # Returns the subset of the pairs given the set of matched ids.
  x = empty.match.result(rke)
  x$match <- subset(rke$pairs, pair.id %in% matched.ids)
  x$status = "OK"
  x$utility = length(matched.ids)
  CHECK_matching(x)
  return(x)
}

get.matching.hospital.utilities <- function(matching, m) {
  hospital.ids <- 1:m  # hospital ids=1,2,3....m
  CHECK_MEMBER(matching$match$hospital, hospital.ids)
  utilities <- sapply(hospital.ids, function(hid) nrow(subset(matching$match, hospital==hid)))
  return(utilities)
}

add.matching <- function(m1, m2) {
  CHECK_matching(m1)
  CHECK_matching(m2)
  CHECK_DISJOINT(m1$match$pair.id, m2$match$pair.id, 
                 msg="Added matchings should have distinct pair ids")
  m = empty.match.result(empty.rke())
  m$match = rbind(m1$match, m2$match)
  m$status = "OK"
  m$utility = m1$utility + m2$utility
  m$information = m1$information + m2$information
  m$matched.cycles = rbind(m1$matched.cycles, m2$matched.cycles)
  return(m)
}

get.matching.ids <- function(matching) {
  CHECK_matching(matching)
  return(matching$match$pair.id)
}

get.matching.utility = function(matching) {
  CHECK_matching(matching)
  return(matching$utility)
}

get.matching.status = function(matching) {
  CHECK_matching(matching)
  return(matching$status)
}

get.matching.3cycles.notAllS <- function(matching) {
  CHECK_matching(matching)
  return(matching$count.3cycles.notallS)
}

empty.match.result <- function(rke) {
  # Empty "matching" object.
  x = subset(rke$pairs, pair.id < 0)
  CHECK_TRUE(nrow(x) == 0, "should be empty")
  ret = list(match=x, status="OK", utility=0, 
             information=empty.matching.information(),
             matched.cycles=matrix(0, nrow=0, ncol=3))
  return(ret)
}

empty.matching.information <- function() {
  # Matches per type.
  # If you edit this, make sure to "match" get.matching.information() as well.
  out = list(info2.OO=0,
             info2.OR=0,
             info2.OS=0,
             info2.OU=0,
             info2.RR=0,
             info2.RS=0,
             info2.RU=0,
             info2.SS=0,
             info2.SU=0,
             info2.UU=0,
             info3.OOO=0,
             info3.OOR=0,
             info3.OOS=0,
             info3.OOU=0,
             info3.ORR=0,
             info3.ORS=0,
             info3.ORU=0,
             info3.OSS=0,
             info3.OSU=0,
             info3.OUU=0,
             info3.RRR=0,
             info3.RRS=0,
             info3.RRU=0,
             info3.RSS=0,
             info3.RSU=0,
             info3.RUU=0,
             info3.SSS=0,
             info3.SSU=0,
             info3.SUU=0,
             info3.UUU=0,
             n2way=0,
             n3way=0)
  ret = as.numeric(out)
  names(ret) <- names(out)
  return(ret)
}

get.matching.information <- function(rke, matched.cycles) {
  # Computes matching information given a matching.
  # Args:
  #   matched.cycles = kx3 matrix with pair ids in the exchange
  # Returns:
  #   An object(LIST) that is the matching information
  #
  # TODO(ptoulis): Add matching information to terminology
  m3 = subset(matched.cycles, type==3)  # matched 3-cycles
  m2 = subset(matched.cycles, type==2)  # swaps 
  m2$type = NULL
  m3$type = NULL
  matched.cycles$type = NULL
  # convention used: infoN.X refers to information on N-way exchange for X info
  # e.g. info3.Sonly = # of S-only 3way exchanges
  #     info2.Ronly = # of R-only swaps (2way exchanges)
  count.type.cycles <- function(cycle.size, pair.types) {
    # Computes |{c}|, c is a cycle that satisfy:
    #  |pairs(c) in TYPES| >= has.at.least
    # i.e.there are "at.least" pairs with a type that belongs to PAIR.TYPES
    CHECK_EQ(cycle.size, length(pair.types))
    m = m2
    if (cycle.size == 3) m = m3
    sum(apply(m, 1, function(xchange.pair.ids) { 
      xchange.pairs = subset(rke$pairs, pair.id %in% xchange.pair.ids); 
      setequal(xchange.pairs$pair.type, pair.types)
    }))
  }
  
  fn <- function(ptypes) count.type.cycles(length(ptypes), pair.types=ptypes)
  out = list(info2.OO=fn(c("O", "O")),
             info2.OR=fn(c("O", "R")),
             info2.OS=fn(c("O", "S")),
             info2.OU=fn(c("O", "U")),
             info2.RR=fn(c("R", "R")),
             info2.RS=fn(c("R", "S")),
             info2.RU=fn(c("R", "U")),
             info2.SS=fn(c("S", "S")),
             info2.SU=fn(c("S", "U")),
             info2.UU=fn(c("U", "U")),
             info3.OOO=fn(c("O", "O", "O")),
             info3.OOR=fn(c("O", "O", "R")),
             info3.OOS=fn(c("O", "O", "S")),
             info3.OOU=fn(c("O", "O", "U")),
             info3.ORR=fn(c("O", "R", "R")),
             info3.ORS=fn(c("O", "R", "S")),
             info3.ORU=fn(c("O", "R", "U")),
             info3.OSS=fn(c("O", "S", "S")),
             info3.OSU=fn(c("O", "S", "U")),
             info3.OUU=fn(c("O", "U", "U")),
             info3.RRR=fn(c("R", "R", "R")),
             info3.RRS=fn(c("R", "R", "S")),
             info3.RRU=fn(c("R", "R", "U")),
             info3.RSS=fn(c("R", "S", "S")),
             info3.RSU=fn(c("R", "S", "U")),
             info3.RUU=fn(c("R", "U", "U")),
             info3.SSS=fn(c("S", "S", "S")),
             info3.SSU=fn(c("S", "S", "U")),
             info3.SUU=fn(c("S", "U", "U")),
             info3.UUU=fn(c("U", "U", "U")),
             n2way=nrow(m2),
             n3way=nrow(m3))
  ret = as.numeric(out)
  names(ret) <- names(out)
  return(ret)
}

gurobi.matched.pairs <- function(gurobi.result, rke, cycles) {
  # Gets the Gurobi output and returns the subset of rke "pairs" object
  # of those that have been matched.
  #
  # Args:
  #   gurobi.result = LIST(status, x=vector of solutions[2cycles, 3cycles])
  # Returns: A <matching> object.
  if(gurobi.result$status=="TIME_LIMIT" | gurobi.result$status=="INF_OR_UNBD") {
    warning("Time limit or infinity.")
    empty.result = empty.match.result(rke)
    empty.result$status = gurobi.result$status
    return(empty.result)
  }
  # The convention is:  x[2-cycles, 3-cycles]
  matched.cycle.ids = which(gurobi.result$x == 1)
  matched.cycles = cycles[matched.cycle.ids, ]
  matched.ids <- c(matched.cycles$pair.id1,
                   matched.cycles$pair.id2, 
                   matched.cycles$pair.id3)
  matched.ids <- setdiff(matched.ids, c(0))  # remove 0's (id3=0 for 2-way xchange)
  not.matched.ids = setdiff(rke.pair.ids(rke), matched.ids)
  CHECK_TRUE(all(!duplicated(matched.ids)), msg="No duplicates in matched ids")
  CHECK_DISJOINT(matched.ids, not.matched.ids, msg="Either matched or not")
 
  # Checks
  CHECK_EQ(length(matched.ids), 2 * nrow(subset(matched.cycles, type == 2)) + 
                                3 * nrow(subset(matched.cycles, type == 3)))
  # Return final result.
  result = subset(rke$pairs, pair.id %in% matched.ids)
  x = empty.match.result(rke)
  x$match = result
  x$status="OK"
  x$utility = length(matched.ids)
  x$information = get.matching.information(rke, matched.cycles)
  x$matched.cycles = matched.cycles
  CHECK_matching(x)
  return (x)
}

max.matching <- function(rke,
                         include.3way=F,
                         ir.constraints=data.frame(),
                         timeLimit=3600,
                         verbose=F) {
  # Performs some form of maximum matching on the specific RKE object.
  #
  # Args:
  #   rke = the RKE object
  #   include.3way = {T, F} whether to use 3-cycles or not
  #   ir.constraints = data.frame with columns (pc, hospital, internal.matches)
  #     Each row has a the #matches of pair "pc" that "hospital" can match internally.
  #
  # Returns:
  #   A MATCHING object. See terminology and above for definitions.
  #   See empty.match.result() for details on the MATCHING object.
  CHECK_rke(rke)
  num.pairs = rke.size(rke)
  num.edges = length(rke.edge.ids(rke))
  # Define Gurobi model 
  # Gurobi defines the problem as: 
  # A * x   <sense>   rhs   ,  sense in {"<=", ">="}
  Cycles = rke.cycles(rke, include.3way=include.3way)
  model.w <- Cycles$type
  if (length(model.w) == 0 | (num.edges == 0)) {
    logthis("Empty RKE", verbose)
    return (empty.match.result(rke))
  }
  model.rhs <- rep(1, num.pairs)
  model.sense <- rep("<=", num.pairs)
  # Build A matrix of model
  # (pairs x cycles) membership matrix. A * x = #cycles/pair (has to be <= 1)
  model.A <- rke.cycles.membership(rke, Cycles)
  get.sub.A <- function(pair.ids) {
    # Get the sub-matrix for only those specific ids.
    if (length(pair.ids) == 0)
      return(matrix(0, nrow=0, ncol=ncol(model.A)))
    ids = as.numeric(rownames(model.A))
    CHECK_MEMBER(pair.ids, ids)
    CHECK_UNIQUE(pair.ids)
    index = match(pair.ids, ids)
    B = matrix(model.A[index,], nrow=length(pair.ids), ncol=ncol(model.A))
    rownames(B) <- pair.ids
    return (B)
  }
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
  
  # Set IR constraints here
  # keep only > 0 constraints
  if (nrow(ir.constraints) > 0) {
    ir.constraints <- subset(ir.constraints, internal.matches > 0)
    if(nrow(ir.constraints) > 0) {
      # 1:0  is a very BAD idea in R.
      for(i in 1:nrow(ir.constraints)) {
        constraint = ir.constraints[i, ]
        # Unload the constraint
        con.pc = constraint$pc
        con.hid = constraint$hospital
        num.matches = constraint$internal.matches
        CHECK_TRUE(!is.null(con.pc) & !is.null(con.hid) & num.matches > 0)
        # For every hospital and pair code
        # 1. Get the pairs subset for hospital=hid, pc=pair code
        sub.pairs = subset(rke$pairs, hospital==con.hid & con.pc==pc)
        pc.hid.pairs <- as.vector(sub.pairs$pair.id)
        CHECK_TRUE(length(pc.hid.pairs) >= num.matches, msg="avail pairs >= matches")
        # 2. Get the pairs x cycles sub-matrix.
        subA = get.sub.A(pc.hid.pairs)
        # 3. Get the (cycles x 1) vector of how many #pairs (hid, pc)
        #   each cycle has (obtained through summing over columns)
        cycle.contrib = as.vector(colSums(subA))
        # Add the constraint
        model$A <- rbind(model$A, cycle.contrib)
        model$sense <- c(model$sense, ">=")
        model$rhs <- c(model$rhs, num.matches)
      }
    }
  }
  
  rownames(model$A) <- 1:nrow(model$A)
  model$A = as.matrix(model$A)
  ##  Seems to be much faster than the old params.
  params.new <- list(OutputFlag=0,                    
                     Cuts=3,
                     Presolve=1,
                     MIPFocus=2,
                     TimeLimit=timeLimit)
  ##  Make sure results > 0
  dimx <- length(model.w)
  model$A <- rbind(model$A, diag(dimx))
  model$rhs <- c(model$rhs, rep(0, dimx))
  model$sense <- c(model$sense, rep(">=", dimx))
  
  gurobi.result <- gurobi(model, params.new)
  # Iron: Just make sure 1.0000   is 1
  gurobi.result$x <- sapply(gurobi.result$x, function(i) {
    if(abs(i) < 1e-6) return(0)
    if(abs(i-1) < 1e-6) return(1)
    return(NA)
  })
  logthis(gurobi.result, verbose)
  match.out = gurobi.matched.pairs(gurobi.result, rke, Cycles)
  edge.index = which(rke$edges$edge.id %in% match.out$matched.edges)
  rke$edges$edge.color[edge.index] <- rep("red", length(match.out$matched.edges))
  #plot.rke(rke)
  if(all(is.element(gurobi.result$x, c(0,1)))) {
    return(match.out)
  } else {
      # Sometimes GUROBI gives results that are not in (0,1)
      # although we have set the variables to be binary (maybe i am missing something)
      # For that reason we call the function again. Empirically, there is no danger of 
     # infinite recursion although this is something that needs to be fixed.
      logwarn("Gurobi unstable output. Saving problematic RKE object AND retrying..")
      save(rke, file="debug/unstable.Rdata")
      print(gurobi.result$x)
      return (max.matching(rke=rke, include.3way=include.3way,
                           ir.constraints=ir.constraints, 
                           timeLimit=timeLimit))
  }
}