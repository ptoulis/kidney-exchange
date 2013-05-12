##   Panos Toulis  ptoulis@fas.harvard.edu
## Contains code implementing different KPD mechanisms. 
##  Currently supported:   rCM,  xCM
## A mechanism receives a combined donor-patient graph and outputs a matching. 
## (note: code returns a vector of utilities (matches/hospital)


## Given the matching, and the input graphs, 
##  calculate the utilities/hospital 
##  Returns    k x 1  vector,     with the matches for every hospital.
get.hospitals.utility <- function(rke.all, m.all) {
    ## find no. of hospitals 
    m = length(unique(rke.all$hospital))
    vals = matrix(0, nrow=m, ncol=1)
    ## get the pair ids in the match.
    matched.ids = m.all$matching$matched.ids
    ## get the hospitals for every pair.
    matched.hospitals =rke.all$hospital[matched.ids]
    
    ## Count matches/hospital. TO-DO(ptoulis): Use "table" here.
    for(hid in 1:m) {
       vals[hid,1] = length(which(matched.hospitals==hid)) 
    }
    return(vals)
}

##  Given an RKE and a type of deviation, 
##  determine which pairs you will *hide*
play.strategy <- function(rke, type) {
  if(type=="t")
    return(c())
  if(type=="c") {
    ## TO-DO(ptoulis): If the max matching is timed-out, this will return empty match
    ## which is equivalent to being truthful. Needs a fix?
    m = max.matching(rke)
    return(m$matching$matched.ids)
  }
  if(type=="r") {
    ## TO-DO(ptoulis): Implement the long-R attack.
    stop("Long-R-attack not implemented.")
  }
  stop("Not Implemented")
}

read.strategy.str = function(strategy.str) {
  strategies = strsplit(strategy.str,split="")[[1]]
  if(length(intersect(strategies, c("t","c", "r")))==0) {
    stop("Strategies should be t,c,r")
  }
  return(strategies)
}
## Implementation of rCM
## deviating: list of hospitals which deviate.
##
##  Return:  kx1   matrix of utilities
rCM <- function(rke.list, strategy.str) {
    m = length(rke.list)
    HospitalUtility = matrix(0, nrow=m, ncol=1)
    
   strategies = read.strategy.str(strategy.str)
    
    for(hid in 1:m) {
      rke.h = rke.list[[hid]]
      matched.internally = play.strategy(rke.h, type=strategies[hid])
      HospitalUtility[hid,1] = length(matched.internally)
      rke.list[[hid]] = remove.pairs(rke.h, matched.internally)
    }
    
    ## 0.  Pool all the reports.
    rke.all = pool.rke(rke.list=rke.list)
    
    ## 1. Simply calculate a maximum-matching (this will shuffle the edges by default)
    m.all =  max.matching(rke.all)
    
    # 2. Compute the utility.
    HospitalUtility = HospitalUtility + get.hospitals.utility(rke.all, m.all)
    
    return(HospitalUtility)
}

## Implementation of  xCM mechanism (Parkes & Toulis, 2013)
##  
##
##
## z = mx1  demand from agents
#  x = supply
g.share = function(z, x) {
  set.J = which(z>0)
  y = rep(0, length(z))
  
  x.r = x 
  while(x.r>0 && length(set.J)>0) {
    if(x.r >= length(set.J)) {
      y[set.J] = y[set.J]+1
      x.r = x.r - length(set.J)
      set.J = setdiff(set.J, which(y>=z))
    } else {
      j.subset = sample(set.J, size=x.r, replace=F)
      y[j.subset] = y[j.subset]+1
      x.r=0
    }
  }
  return(y)
}
xCM <- function(rke.list, strategy.str) {
   # total no. of hospitals
    m = length(rke.list)
    strategies = read.strategy.str(strategy.str)
    ## Utility = K x 1 {  total matched }
    HospitalUtility = matrix(0, nrow=m, ncol=1)
    
    ##  Initialize variables
    ## rke.list.S = S-subgraphs of hospitals
    rke.list.S = list()
    rke.list.R = list()
    ##  Demand for R-pairs
    z.AB = c()
    z.BA = c()
    ## IR constraints per hospital
    IR.constraints.S = list()
    IR.constraints.R = list()
    
    ## Pair code (useful when setting constraints)
    pc.AB = pair.code(list(donor="A", patient="B"))
    pc.BA = pair.code(list(donor="B", patient="A"))
    pc.R = c(pc.AB, pc.BA)
    ## Iterate over all hospitals.
    for(hid in 1:m) {
      ## Constraints.
        irs = rep(0, length(Pair.Codes))
        irr = rep(0, length(Pair.Codes))
        
        # 0. Play the strategy first.
        ##  update rke.h   and the rke list 
        rke.h = rke.list[[hid]]
        matched.internally = play.strategy(rke.h, type=strategies[hid])
        rke.h = remove.pairs(rke.h, matched.internally)
        rke.list[[hid]]  = rke.h
        ## Matched internally, so add the utility now.
        HospitalUtility[hid,1] = length(matched.internally)
        
        ## xCM mechanism starts
        #  1. Compute max matching internally in S
        S.subgraph = get.subgraph(rke.h, type="S")
        mS = max.matching(S.subgraph)
        rke.list.S[[hid]] = S.subgraph
        matches.tab = table(rke.h$pc[mS$matching$matched.ids])
        pcs.matched = as.numeric(names(matches.tab))
        no.matches = as.numeric(matches.tab)
        # S-subgraph constraints
        irs[pcs.matched] = no.matches
        
        
        # 2. Matchings in R 
        R.subgraph = get.subgraph(rke.h, type="R")
        mR = max.matching(R.subgraph)
        rke.list.R[[hid]] = R.subgraph
        # Matched A-B, B-A pairs
        N.AB = length(intersect(filter.pairs.by.donor.patient(rke.h, dtype="A",ptype="B"),
                                mR$matching$matched.ids))
        N.BA = length(intersect(filter.pairs.by.donor.patient(rke.h, dtype="B",ptype="A"),
                                mR$matching$matched.ids))
        ## Because we are working on the R-subgraph, these numbers should be the same.
        if(N.AB != N.BA) 
          stop("AB and BA matches should be equal! ")
        ## # unmatched A-B, B-A pairs
        z.AB[hid] = length(intersect(filter.pairs.by.donor.patient(rke.h, dtype="A",ptype="B"),
                                     mR$matching$not.matched.ids))
        z.BA[hid] = length(intersect(filter.pairs.by.donor.patient(rke.h, dtype="B",ptype="A"),
                                     mR$matching$not.matched.ids))
        
        ## add R constraints.
        irr[pc.R] = c(N.AB, N.BA)
        
        
        IR.constraints.S[[hid]] = irs
        IR.constraints.R[[hid]] = irr
    }
    ###  Done with per-hospital
    ###  Final stages of xCM
    x.AB = sum(z.AB)
    x.BA = sum(z.BA)
    y.AB = rep(0, m)
    if(x.AB>=x.BA)
      y.AB = g.share(z.AB, x.BA)
    
    y.BA = rep(0,m)
    if(x.AB<x.BA)
      y.BA = g.share(z.BA, x.AB)
    
    ###  1.  Match S internally
    rke.all = pool.rke(rke.list)
    all.but.s = setdiff(all.edges(rke.all), filter.edges.by.type(rke.all, "S","S"))
    match.s = max.matching(rke.all, IR.constraints=IR.constraints.S,
                           remove.edges = all.but.s)
    
    ## 2.   Match R internally
    all.but.r = setdiff(all.edges(rke.all), filter.edges.by.type(rke.all, "R","R"))
    Kq = c()
    match.r = list()
    q = 0
    while(length(Kq)==0) {
      ir.constraints = list()
      for(h in 1:m)
        ir.constraints[[h]][pc.AB] = IR.constraints.R[[h]][pc.AB]+max(0, y.AB-q)
      for(h in 1:m)
        ir.constraints[[h]][pc.BA] = IR.constraints.R[[h]][pc.BA]+max(0, y.BA-q)
      
      ## Do the matching.
      match.r = max.matching(rke.all, IR.constraints = ir.constraints,
                             remove.edges = all.but.r)
      Kq = match.r$matching$matched.ids
      q = q + 1
    }
    
    #  3. Almost regular matching to the remainder
    remainder = remove.pairs(rke.all, union(match.r$matching$matched.ids, 
                                            match.s$matching$matched.ids))
    match.od = max.matching(remainder, regular.matching=T)
    
    HospitalUtility = HospitalUtility + get.hospitals.utility(rke.all, match.s)
                                      + get.hospitals.utility(rke.all, match.r)
                                      + get.hospitals.utility(remainder, match.od)

    return(HospitalUtility)
}