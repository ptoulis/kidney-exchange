##   Panos Toulis  ptoulis@fas.harvard.edu
## Contains code implementing different KPD mechanisms. 
##  Currently supported:   rCM,  xCM
## A mechanism receives a combined donor-patient graph and outputs a matching. 
## (note: code returns a vector of utilities (matches/hospital)
source("rke.R")
source("matching.R")

Sample.Setup <- function(m, n, strategy.str, uniform.pra=T) {
  
  rke.list = rrke.many(m=m, n=n, uniform.pra= uniform.pra)
  rke.all = pool.rke(rke.list)
  
  x = init.mechanism(rke.list, strategy.str)
  setup = list()
  
  reported.rke.list = list()
  ## hidden.pairs = ids in terms of rke.all  that were hidden
  hidden.pairs = c()
  #  1b.  Populate the reported and hidden graphs.
  for(h in 1:m) {
    x.h = x[[h]]
    reported.rke.list[[h]] = remove.pairs(rke.list[[h]], 
                                          x.h$hide)
    hidden.h = x.h$hide
    if(length(hidden.h)>0)
        hidden.pairs = c(hidden.pairs, 
                           sapply(hidden.h, function(i) 
                                            rke.all$map.fwd(i, h)))
  }
  ##  Save the "setup" object
  ## This is important in order to compare mechanisms
  ## on the same sets of hospitals and pooled graphs.
  setup$reported.rke.list = reported.rke.list
  setup$reported.rke.all = remove.pairs(rke.all, hidden.pairs)
  setup$rke.list = rke.list
  setup$rke.all = rke.all
  return(setup)
  
}
##  Runs a mechanism
Run.Mechanism = function(setup, mech) {
  
  reported.rke.list = setup$reported.rke.list
  reported.rke.all = setup$reported.rke.all
  rke.list = setup$rke.list
  rke.all = setup$rke.all
  
  ## 3. Run the mechanism
  mech.out.ids = do.call(mech, args=list(rke.list=reported.rke.list, 
                                         rke.all = reported.rke.all) )

  ## 4. Compute utility from mechanism
  Util = get.hospitals.utility(reported.rke.all, mech.out.ids)

  ## 5.  Utility from final internal matches. 
  ##     
  for(h in 1: length(rke.list)) {
    rke.h = rke.list[[h]]
    matched.ids.in.reported.all = intersect( get.hospital.pairs(reported.rke.all, h), 
                                  mech.out.ids )
    matched.ids.in.all = reported.rke.all$map.back(matched.ids.in.reported.all)
    matched.ids.in.hosp = rke.all$map.back(matched.ids.in.all)
    
    rke.remainder = remove.pairs(rke.h, matched.ids.in.hosp)
    m.h = max.matching(rke.remainder)
    Util[h] = Util[h] + m.h$matching$utility

  }
  return(Util)
}


## Given the matching, and the input graphs, 
##  calculate the utilities/hospital 
##  Returns    k x 1  vector,     with the matches for every hospital
get.hospitals.utility <- function(rke.all, matched.ids) {
    ## find no. of hospitals 
    m = length(unique(rke.all$hospital))
    vals = matrix(0, nrow=m, ncol=1)
   
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
# Returns:    list(  hide= c(..),   report=c(...) )
play.strategy <- function(rke, type) {
  ret = list(report=c(), hide=c())
  
  if(type=="t") {
    ret$hide = c()
  } else if(type=="c") {
    ## TO-DO(ptoulis): If the max matching is timed-out, this will return empty match
    ## which is equivalent to being truthful. Needs a fix?
    m = max.matching(rke)
    ret$hide = m$matching$matched.ids
  }  else if(type=="r") {
    pairs.AB = filter.pairs.by.donor.patient(rke, "A","B")
    pairs.BA = filter.pairs.by.donor.patient(rke, "B","A")
    pairs.O = filter.pairs.by.type(rke,"O")
    ## want to hide the "short" side of R
    hide.R= pairs.AB
    if(length(pairs.BA) < length(pairs.AB))
      hide.R = pairs.BA
    
     ret$hide = c(hide.R)
  }
  ret$report = setdiff(1:get.size(rke), ret$hide)
  return(ret)
}

read.strategy.str = function(strategy.str) {
  strategies = strsplit(strategy.str,split="")[[1]]
  if(length(intersect(strategies, c("t","c", "r")))==0) {
    stop("Strategies should be t,c,r")
  }
  return(strategies)
}

## Initializes every mechanism + plays strategy
## Returns:  list(  hid => play.strategy )
init.mechanism = function(rke.list, strategy.str) {

  m = length(rke.list)  
  strategies = read.strategy.str(strategy.str)
  strategy.list = list()
  
  for(hid in 1:m) {
    rke.h = rke.list[[hid]]
    strategy.list[[hid]] = play.strategy(rke.h, type=strategies[hid])
  }
  return( strategy.list )
}

## Implementation of rCM
## Returns:   c(...)   = matched ids
rCM <- function(rke.list, rke.all) {

  ## 1. Simply calculate a maximum-matching (this will shuffle the edges by default)
  m.all =  max.matching(rke.all)
  
  matched.all.ids = m.all$matching$matched.ids
  # 2. Compute the utility.
  return(matched.all.ids)
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

## TO-DO(ptoulis): Don't really care much about correctness
## as long as it has good properties, we can treat it as a black box.
xCM <- function(rke.list, rke.all) {
  warning("xCM() has no unit-test")
  
  m = length(rke.list)
  matched.all.ids = c()
  
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
                     match.s$matching$matched.ids, "R and S matched ids")
  
  redges = match.r$matching$matched.edges
  sedges = match.s$matching$matched.edges
  
  TEST.SETS.DISJOINT(redges, sedges, "S and R edges")
  ###  Update the matching output
  matched.all.ids = union(match.r$matching$matched.ids, 
                          match.s$matching$matched.ids)
  matched.edges.already = union( redges, sedges )
  
  ## Match OD's individually.
  for(hid in 1:m) {
    Gh = get.hospital.pairs(rke.all, hid)
    Gh.remainder = setdiff(Gh, intersect(Gh, matched.all.ids))
    
    match.Gh = max.matching(rke.all, regular.matching=T, 
                            remove.edges= get.external.edges(rke.all, Gh.remainder))
    
    ## Make sure OD ids are not R or S pairs
    TEST.SETS.DISJOINT(match.Gh$matching$matched.ids, 
                       matched.all.ids, str="OD and matched_already pairs")
    TEST.SETS.DISJOINT(match.Gh$matching$matched.edges, matched.edges.already,
                       str="Matched edges already")
    
    matched.Gh.ids = match.Gh$matching$matched.ids
    matched.all.ids = c(matched.all.ids,  matched.Gh.ids)
  }
  return(matched.all.ids)
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
                      QS) {
  
  #warning("UD lottery does not have a unit test.")
  Qh = QS$Q
  Sh = QS$S
  
  # The underdemanded PC code
  ud.code = pair.code(pair.ud)
  
  # X-Y pairs in rke
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
    
    Bh.XY[[hid]] = intersect(rke.XY.pairs, get.hospital.pairs(rke, hid) )
    ## Sanity check
    TEST.SUBSET(Sh.XY[[hid]], Bh.XY[[hid]])
  }
  
  get.sum = function() sum(sapply(Hn, function(hid) length(Sh.XY[[hid]])))
  
  plateau.count = 0
  
  while(get.sum() < theta && 
          plateau.count<10 && 
          length(set.J)>0) {
    h.sample = uniform.sample(set.J)
    set.J = set.J[-which(set.J==h.sample)[1]]
    
    TEST.SUBSET(Sh.XY[[h.sample]], Bh.XY[[h.sample]] )
    
    avail = setdiff(Bh.XY[[h.sample]], Sh.XY[[h.sample]])
    if(length(avail)>0) {
      xy.pair = uniform.sample(avail)
      Sh.XY[[h.sample]] = c(Sh.XY[[h.sample]], xy.pair)
      plateau.count=0
    } else {
      plateau.count = plateau.count +1 
    }
  }
  
  ret.pairs = c()
  # return the union.
  for(h in Hn) {
    ret.pairs = c(ret.pairs, Sh.XY[[h]])
  }
  return(ret.pairs)
}

Bonus.QS = function(rke.all)  {
  
  #  Before we start. Compute Qh and Sh sets internally
  ## Q[hid][X-Y] = how many #X-Y  in hospital hid 
  ## S[hid][X-Y] = { ids } of type X-Y which are included in a max regular matching.
  Qh = list()
  Sh = list()
  all.hospitals = sort( unique(rke.all$hospital) )
  ud.pcs = pair.codes.per.type("U")
  
  for(h in all.hospitals) {
    Qh[[h]] = rep(0, length(Pair.Codes))
    Sh[[h]] = list()
    h.pairs= get.hospital.pairs(rke.all, h)
    h.pcs = rke.all$pc[h.pairs]
    # 2 . Compute a regular matching on the specific hospital (internal)
    not.h.edges = get.external.edges(rke.all, pair.ids=h.pairs)
    
    m = max.matching(rke.all, regular.matching=T, 
                     remove.edges= not.h.edges )
    ## Iterate over all UD codes
    for(i in ud.pcs) {
      hpairs.of.type = h.pairs[ which(h.pcs==i) ]
      Qh[[h]][i] = length(hpairs.of.type)
      Sh[[h]][[i]] = intersect( m$matching$matched.ids,  hpairs.of.type)
    }
  }
  ##  TO-DO(ptoulis): How do you know if you computed Qh, Sh correctly???
  ret = list(Q=Qh, S=Sh)
}
## Bonus mechanism. Ashlagi & Roth (2013)
Bonus = function(rke.list, rke.all) {
  
  matched.all.ids = c()
  ## 0. Initialize mechanism
  m = length(rke.list)
  
  ## Pair code (useful when setting constraints)
  pc.AB = pair.code(list(donor="A", patient="B"))
  pc.BA = pair.code(list(donor="B", patient="A"))
  pc.R = c(pc.AB, pc.BA)
  
  IR.constraints = compute.ir.constraints(rke.list, types=c("R"))
  
  ## 1. Match S pairs
  all.but.s = get.external.edges( rke.all, filter.pairs.by.type(rke.all, "S") )
  match.s = max.matching(rke.all, remove.edges = all.but.s)
  
  ## 2. Match R pairs
  all.but.r = get.external.edges( rke.all, filter.pairs.by.type(rke.all, "R") )
  match.r = max.matching(rke.all, IR.constraints=IR.constraints$R, 
                         remove.edges = all.but.r)
  
  #ORr = length(which ( get.pair.types(rke.all, match.r$matching$matched.ids) == "O" ) )
  #ORs = length(which ( get.pair.types(rke.all, match.r$matching$matched.ids) == "O" ) )
  
  TEST.SUBSET( get.pair.types(rke.all, match.r$matching$matched.ids), c("R"), str=" only R")
  ## Check if matched.already containts *only* R-S pairs.
  TEST.SUBSET(get.pair.types(rke.all, match.s$matching$matched.ids), c("S"), str="only S pairs")
  
  matched.ids.R = match.r$matching$matched.ids
  matched.ids.S = match.s$matching$matched.ids
  
  TEST.SETS.DISJOINT(matched.ids.R, matched.ids.S)
  
  matched.all.ids = c(matched.ids.R, matched.ids.S)
  
  #  3. Match OD/UD pairs.
  k = as.integer(m/2)
  H.sets = list()
  H.sets[[1]] = sample(1:m, size=k, replace=F)
  H.sets[[2]] = setdiff(1:m, H.sets[[1]])
  
  ## Underdemanded PC codes
  ud.pc = pair.codes.per.type("U")
  # Compute Q-S, Q[hid][X-Y] = how many #X-Y in hospital hid
  # S[hid][X-Y] = {}  ids of X-Y in hospital hid 
  QS.obj = Bonus.QS(rke.all)
  
  
  ## For all under-demanded pairs  X-Y
  for(i in ud.pc) {
    
    pair.ud = pair.code.to.pair(i)
    # Find the reciprocal OD pair  Y-X
    pair.od = list(donor=pair.ud$patient, patient=pair.ud$donor)
    ## Note: follow the notation in Ashlagi&Roth about Bonus
    ## Find all pairs Y-X
    all.YX = filter.pairs.by.donor.patient(rke.all, pair.od$donor, pair.od$patient)
    all.XY = filter.pairs.by.donor.patient(rke.all, pair.ud$donor, pair.ud$patient)
    
    hospitals.YX = rke.all$hospital[all.YX]
    ## Iterate through hospital groups
    for(j in 1:2) {
      
      ## Pairs that belong to the other group
      BHother.pairs =  which(rke.all$hospital %in% H.sets[[3-j]])
      BHj.pairs     =  which(rke.all$hospital %in% H.sets[[j]])

      ## Pairs Y-X in this set of hospitals
      tau.BHother.YX =  intersect(all.YX,  BHother.pairs)
      tau.BHj.XY =  intersect(all.XY,  BHj.pairs)
      
      ## no. of such pairs
      theta.j.YX = length( tau.BHother.YX )
      
      ## Run the under-demanded lottery  (contains pair ids)
      S.XY = ud.lottery(rke.all, pair.ud, H.sets[[j]], theta.j.YX, QS.obj)
      
      # Test if disjoint (former=UD latter=OD so should not intersect)
      TEST.SETS.DISJOINT(S.XY, tau.BHother.YX, str="S.XY with BH.XY")
      TEST.SUBSET(S.XY, tau.BHj.XY, str="S.XY subset")
      ## the only ids to be considered in the matching
      xyyx.ids = union(S.XY, tau.BHother.YX)
      exclude.edges = get.external.edges(rke.all, pair.ids= xyyx.ids)
      
      ##  Here you match   X-Y pairs from Hj   with   Y-X pairs from Hother
      Mj = max.matching(rke.all, regular.matching=F, 
                           remove.edges= exclude.edges)
      ## Test if newly matched have already been matched.
      TEST.SETS.DISJOINT(Mj$matching$matched.ids, matched.all.ids, str="Matched now vs. already")
      
      matched.all.ids = c(matched.all.ids, Mj$matching$matched.ids)

    }
  }
  return(matched.all.ids)
}

