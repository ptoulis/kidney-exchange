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

## Initializes every mechanism. Do not change.
init.mechanism = function(rke.list, strategy.str) {

  m = length(rke.list)
  HospitalUtility = matrix(0, nrow=m, ncol=1)
  
  strategies = read.strategy.str(strategy.str)
  
  for(hid in 1:m) {
    rke.h = rke.list[[hid]]
    matched.internally = play.strategy(rke.h, type=strategies[hid])
    HospitalUtility[hid,1] = length(matched.internally)
    rke.list[[hid]] = remove.pairs(rke.h, matched.internally)
  }
  return( list(rke.list=rke.list, util=HospitalUtility)   )
}
## Implementation of rCM
## deviating: list of hospitals which deviate.
##
##  Return:  kx1   matrix of utilities
rCM <- function(rke.list, strategy.str) {
    warning("rCM does not have unit test")
    x = init.mechanism(rke.list, strategy.str)
    rke.list = x$rke.list
    HospitalUtility = x$util
    
    ## 0.  Pool all the reports.
    rke.all = pool.rke(rke.list=rke.list)
    
    ## 1. Simply calculate a maximum-matching (this will shuffle the edges by default)
    m.all =  max.matching(rke.all)
    
    # 2. Compute the utility.
    HospitalUtility = HospitalUtility + get.hospitals.utility(rke.all, m.all)
    
    return(HospitalUtility)
}


## Compute IR constraints.
compute.ir.constraints = function(rke.list, types=c()) {
  m = length(rke.list)
  if(length(types)==0 || m==0 )
    return(list() )
  ##  Initialize variables
  ##  Demand for R-pairs
  z.AB = rep(0, m)
  z.BA = rep(0, m)

  ## IR constraints per hospital
  IR.constraints = list(S=list(), R=list())
  
  ## Pair code (useful when setting constraints)
  pc.AB = pair.code(list(donor="A", patient="B"))
  pc.BA = pair.code(list(donor="B", patient="A"))
  pc.R = c(pc.AB, pc.BA)
  ## Iterate over all hospitals.
  for(hid in 1:m) {
    ## Constraints.
    irs = rep(0, length(Pair.Codes))
    irr = rep(0, length(Pair.Codes))
    
    rke.h = rke.list[[hid]]
    S.subgraph = get.subgraph(rke.h, type="S")
    R.subgraph = get.subgraph(rke.h, type="R")
    
    ## xCM mechanism starts
    #  1. Compute max matching internally in S
    if("S" %in% types) {
      mS = max.matching(S.subgraph)
      if(length(mS$matching$matched.ids)>0) {
        matches.tab = table(S.subgraph$pc[mS$matching$matched.ids])
        pcs.matched = as.numeric(names(matches.tab))
        no.matches = as.numeric(matches.tab)
        irs[pcs.matched] = no.matches
      } 
      # done S-subgraph constraints
    }
    if("R" %in% types) {
      # 2. Matchings in R 
      mR = max.matching(R.subgraph)
      # Matched A-B, B-A pairs
      N.AB = length(intersect(filter.pairs.by.donor.patient(R.subgraph, dtype="A",ptype="B"),
                              mR$matching$matched.ids))
      N.BA = length(intersect(filter.pairs.by.donor.patient(R.subgraph, dtype="B",ptype="A"),
                              mR$matching$matched.ids))
      ## Because we are working on the R-subgraph, these numbers should be the same.
      if(N.AB != N.BA) 
        stop("AB and BA matches should be equal! ")
      if(N.AB + N.BA != length(mR$matching$matched.ids)) 
        stop("Something wrong with the subgraph. Not only A-B, B-A pairs considered.!")
      ## # unmatched A-B, B-A pairs
      z.AB[hid] = length(intersect(filter.pairs.by.donor.patient(R.subgraph, dtype="A",ptype="B"),
                                   mR$matching$not.matched.ids))
      z.BA[hid] = length(intersect(filter.pairs.by.donor.patient(R.subgraph, dtype="B",ptype="A"),
                                   mR$matching$not.matched.ids))
      
      ## add R constraints.
      irr[pc.R] = c(N.AB, N.BA)
    }
    
    ## do some cleaning
    rm(list=c("R.subgraph", "S.subgraph", "rke.h"))
    
    IR.constraints$S[[hid]] = irs
    IR.constraints$R[[hid]] = irr
  }## for every hospital
  IR.constraints$z.ab = z.AB
  IR.constraints$z.ba = z.BA
  return(IR.constraints)
}
## Implementation of  xCM mechanism (Parkes & Toulis, 2013)
##  
##
##
## z = mx1  demand from agents
#  x = supply
# TO-DO(ptoulis): Unit tests for this one?
g.share = function(z, x) {
  warning("g.share has no unit-test")
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
    warning("xCM() has no unit-test")
    # 0. Initialize mechanism
    m =length(rke.list)
    x = init.mechanism(rke.list, strategy.str)
    HospitalUtility = x$util
    rke.list = x$rke.list
    
    rke.all = pool.rke(rke.list)
  
    ##  1. Compute IR constraints
    IR.constraints = compute.ir.constraints(rke.list, types=c("S", "R"))
    
    
    z.AB = IR.constraints$z.ab
    z.BA = IR.constraints$z.ba
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
    
    ##  Ready to run xCM now
    ###  2.  Match S internally
    all.but.s = filter.out.edges.by.type(rke.all, "S","S")
    match.s = max.matching(rke.all, IR.constraints=IR.constraints$S,
                           remove.edges = all.but.s)
    
    ## 3.   Match R internally
    ## TO-DO(ptoulis): Slow for some reason
    all.but.r = filter.out.edges.by.type(rke.all, "R","R")
    Kq = c()
    match.r = list()
    q = 0    
    ## Pair code (useful when setting constraints)
    pc.AB = pair.code(list(donor="A", patient="B"))
    pc.BA = pair.code(list(donor="B", patient="A"))
    pc.R = c(pc.AB, pc.BA)
    ## TO-DO(ptoulis): very slow!!
    while(length(Kq)==0) {
      ir.constraints = list()
      for(h in 1:m)
        ir.constraints[[h]] = rep(0, length(Pair.Codes))
      
      for(h in 1:m) {
        ir.constraints[[h]][pc.AB] = IR.constraints$R[[h]][pc.AB]+max(0, y.AB[h]-q)
        ir.constraints[[h]][pc.BA] = IR.constraints$R[[h]][pc.BA]+max(0, y.BA[h]-q)
      }
      
      ## Do the matching.
      match.r = max.matching(rke.all, IR.constraints = ir.constraints,
                             remove.edges = all.but.r)
      Kq = match.r$matching$matched.ids
      q = q + 1
    }
    print(sprintf("Final q* = %d", q))
    ## remove some stuff that are not needed anymore
    rm(IR.constraints)
    
    
    #  4. Almost regular matching to the remainder
    ## Notice that match.r, match.s are all on rke.all so that the 
    ## ids refer to the same original ids in rke.all
    TEST.SETS.DISJOINT(match.r$matching$matched.ids, 
                        match.s$matching$matched.ids)
    
    redges = match.r$matching$matched.edges
    sedges = match.s$matching$matched.edges
    
    TEST.SETS.DISJOINT(redges, sedges)
    
    matched.already = union(match.r$matching$matched.ids, 
                                match.s$matching$matched.ids)
    matched.edges.already = union( redges, sedges )
    
    remainder = remove.pairs(rke.all, matched.already)
    match.od = max.matching(rke.all, regular.matching=T, 
                            remove.edges= get.incident.edges(rke.all, matched.already))
    ## Make sure OD ids are not R or S pairs
    TEST.SETS.DISJOINT(match.od$matching$matched.ids, 
                       matched.already)
    TEST.SETS.DISJOINT(match.od$matching$matched.edges, matched.edges.already)
    
    Us = get.hospitals.utility(rke.all, match.s)
    Ur = get.hospitals.utility(rke.all, match.r)
    Uo = get.hospitals.utility(rke.all, match.od)

    HospitalUtility = HospitalUtility +Ur+ Us+ Uo

    return(HospitalUtility)
}


## The underdemanded lottery as defined in Ashlagi & Roth (2013)
##  Return set of pair ids
## rke = entire RKE object
## pair.ud = UD pair (donor, patient)
# Hn =   [1.4.5....]   vector of hospital ids
# theta = supply
# Qh, Sh  = lists 
## @RETURN: (id1, id2, ...)    list of pair ids which are the union of Sh(X-Y)
ud.lottery = function(rke, 
                      pair.ud, 
                      Hn, 
                      theta,
                      Qh, Sh) {
  warning("UD lottery does not have a unit test.")
  # The underdemanded PC code
  ud.code = pair.code(pair.ud)
  
  rke.XY.pairs = filter.pairs.by.donor.patient(rke, pair.ud$donor, pair.ud$patient)
  ##  Pre-compute the X-Y pairs for every hospital h
  Bh.XY = list()
  ## For every hospital, holds the set of X-Y pairs.
  Sh.XY = list()
  ## Contains # balls /hospital     [1,1,1,3,3,3,3,3,5,5,5, ...]   NOTE: this is index by a hospital in Hn
  set.J = c()
  ## Initialize the Sh.XY  to the internal max matches.
  for(hid in Hn) {
    Sh.XY[[hid]] = Sh[[hid]][[ud.code]]
    set.J = c(set.J, rep(hid, Qh[[hid]][ud.code]))
    
    Bh.XY[[hid]] = intersect(rke.XY.pairs, which(rke$hospital==hid) )
    ## Sanity check
    if(!is.subset(Bh.XY[[hid]], Sh.XY[[hid]]))
      stop("Sh.XY should be subset. Error in UD lottery")
  }

  get.sum = function() sum(sapply(Hn, function(hid) length(Sh.XY[[hid]])))
  
  trials = 0
  
  while(get.sum() < theta && trials<20 && length(set.J)>0) {
    trials= trials+1
    h.sample = uniform.sample(set.J)
    set.J = set.J[-which(set.J==h.sample)[1]]
    if(! is.subset(Bh.XY[[h.sample]], Sh.XY[[h.sample]]))
      stop("Sh should be a subset of Bh! ")
    
    avail = setdiff(Bh.XY[[h.sample]], Sh.XY[[h.sample]])
    if(length(avail)>0) {
      xy.pair = uniform.sample(avail)
      Sh.XY[[h.sample]] = c(Sh.XY[[h.sample]], xy.pair)
    }
  }
  ret.pairs = c()
  # return the union.
  for(h in Hn) {
    ret.pairs = c(ret.pairs, Sh.XY[[h]])
  }
  return(ret.pairs)
}


## Bonus mechanism. Ashlagi & Roth (2013)
Bonus = function(rke.list, strategy.str) {
  
  ## 0. Initialize mechanism
  m = length(rke.list)
  x = init.mechanism(rke.list, strategy.str)
  HospitalUtility = x$util
  rke.list = x$rke.list
  rke.all = pool.rke(rke.list)
  
  ## Pair code (useful when setting constraints)
  pc.AB = pair.code(list(donor="A", patient="B"))
  pc.BA = pair.code(list(donor="B", patient="A"))
  pc.R = c(pc.AB, pc.BA)
  
  IR.constraints = compute.ir.constraints(rke.list, types=c("R"))
  
  ## 1. Match S pairs
  all.but.s = filter.out.edges.by.type(rke.all, "S","S")
  match.s = max.matching(rke.all, remove.edges = all.but.s)
    
  ## 2. Match R pairs
  all.but.r = filter.out.edges.by.type(rke.all, "R","R")
  match.r = max.matching(rke.all, IR.constraints=IR.constraints$R, 
                         remove.edges = all.but.r)
  
  matched.ids.R = match.r$matching$matched.ids
  matched.ids.S = match.s$matching$matched.ids
  
  TEST.SETS.DISJOINT(matched.ids.R, matched.ids.S)
  
  matched.already = c(matched.ids.R, matched.ids.S)
 
  
  Us = get.hospitals.utility(rke.all, match.s)
  Ur = get.hospitals.utility(rke.all, match.r)
  
  ## Update utilities.
  HospitalUtily = HospitalUtility + Ur + Us
  ##   Standard up to here.
  
  #  3. Match OD/UD pairs.
  k = as.integer(m/2)
  H.sets = list()
  H.sets[[1]] = sample(1:m, size=k, replace=F)
  H.sets[[2]] = setdiff(1:m, H.sets[[1]])

  ## Underdemanded PC codes
  ud.pc = which( sapply(Pair.Codes, function(i) pair.type(pair.code.to.pair(i)))=="U")
  
  #  Before we start. Compute Qh and Sh sets internally
  ## Q[hid][X-Y] = how many #X-Y  in hospital hid 
  ## S[hid][X-Y] = { ids } of type X-Y which are included in a max regular matching.
  Qh = list()
  Sh = list()
  for(h in 1:m) {
    Qh[[h]] = rep(0, length(Pair.Codes))
    Sh[[h]] = list()
    h.pairs= which(rke.all$hospital == h)
    h.pcs = rke.all$pc[h.pairs]
    # 2 . Compute a regular matching on the specific hospital (internal)
    edges.score = sapply(rke.edges(rke.all), function(e) 
                  sum(get.incident.nodes(rke.all, c(e)) %in% h.pairs))
    not.h.edges = which(edges.score<2)
    
    m = max.matching(rke.all, regular.matching=T, 
                     remove.edges= not.h.edges )
    for(i in ud.pc) {
      hpairs.of.type = h.pairs[ which(h.pcs==i) ]
      Qh[[h]][i] = length(hpairs.of.type)
      Sh[[h]][[i]] = intersect( m$matching$matched.ids,  hpairs.of.type)
    }
  }
  ##  TO-DO(ptoulis): How do you know if you computed Qh, Sh correctly???
  
  ## For all under-demanded pairs  X-Y
  for(i in ud.pc) {
    
    pair.ud = pair.code.to.pair(i)
    # Find the reciprocal OD pair  Y-X
    pair.od = list(donor=pair.ud$patient, patient=pair.ud$donor)
    ## Note: follow the notation in Ashlagi&Roth about Bonus
    ## Find all pairs Y-X
    set.YX = filter.pairs.by.donor.patient(rke.all, pair.od$donor, pair.od$patient)
    hospitals.YX = rke.all$hospital[set.YX]
    ## Iterate through hospital groups
    for(j in 1:2) {
     
      ## Pairs that belong to the other group
      BH =  which(rke.all$hospital %in% H.sets[[3-j]])
      ## Pairs Y-X in this set of hospitals
      tau.BH.YX =  set.YX[ which(hospitals.YX %in% BH) ] 
      ## no. of such pairs
      theta.j.YX = length( tau.BH.YX )
      
      ## Run the under-demanded lottery  (contains pair ids)
      S.XY = ud.lottery(rke.all, pair.ud, H.sets[[j]], theta.j.YX,
                        Qh, Sh)
      
      # Test if disjoint (former=UD latter=OD so should not intersect)
      TEST.SETS.DISJOINT(S.XY, tau.BH.YX)
      
      inc.ids = union(S.XY, tau.BH.YX)
      exclude.edges = get.incident.edges(rke.all, matched.already)
      set.remove.edges =  union(exclude.edges, get.nonincident.edges(rke.all, inc.ids))
      Mj.XY = max.matching(rke.all, regular.matching=F, 
                           remove.edges= set.remove.edges)
      ## Test if newly matched have already been matched.
      TEST.SETS.DISJOINT(Mj.XY$matching$matched.ids, matched.already)
      
      matched.already = c(matched.already, Mj.XY$matching$matched.ids)
      
      HospitalUtility = HospitalUtility + get.hospitals.utility(rke.all, Mj.XY)
    }
  }
  
  return(HospitalUtility)
}

