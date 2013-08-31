# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Contains code implementing different KPD mechanisms. 
##  Currently supported:   rCM,  xCM
## A mechanism receives a combined donor-patient graph and outputs a matching. 
## (note: code returns a vector of utilities (matches/hospital)
source("rke.R")
source("matching.R")

kpd.create <- function(rke.pool, strategy.str, verbose=F) {
  # Creates a KPD from the specific pool.
  # Proceeds as follows:
  # (1) Plays the strategy for every hospital
  # (2) Creates a sub-rke object by removing all pairs that hospitals hide
  # (3) Creates a KPD object from the reported RKE pool and real RKE pool (input)
  #
  # Args:
  #   rke.pool : The real rke pool object (what the hospitall really own)
  #   strategy.str : Char vector of strategies (one of "t,b,c,r")
  #                  Has to be same #length as #hospitals
  # Returns:
  #   The KPD object from the reported and real RKE pools.
  CHECK_rke.pool(rke.pool)
  rke.list = rke.pool$rke.list
  rke.all = rke.pool$rke.all
  m = length(rke.list)
  CHECK_EQ(nchar(strategy.str), m, msg="Correct #of strategies")
  logthis("Playing strategies", verbose)
  hospital.action = play.strategies(rke.list, strategy.str)
  reported.rke.list = list()
  ## hidden.pairs = ids in terms of rke.all  that were hidden
  all.hidden.pairs = c()
  #  1b.  Populate the reported and hidden graphs.
  for(hid in 1:m) {
    hospitalAction = hospital.action[[hid]]
    logthis(sprintf("Hospital %d is hiding %d pairs", hid, length(hospitalAction$hide)),
            verbose)
    reported.rke.list[[hid]] = rke.remove.pairs(rke.list[[hid]],
                                                hospitalAction$hide)
    all.hidden.pairs = c(all.hidden.pairs, hospitalAction$hide)
  }
  ##  Save the "kpd" object
  ## This is important in order to compare mechanisms
  ## on the same sets of hospitals and pooled graphs.
  kpd = list()
  kpd$reported.pool = list(rke.list=reported.rke.list,
                           rke.all=rke.remove.pairs(rke.all, all.hidden.pairs))
  kpd$real.pool = list(rke.list=rke.list, rke.all = rke.all)
  CHECK_kpd(kpd)
  return(kpd)
}

play.strategy <- function(rke, strategy.str, include.3way=F) {
  # For the RKE, play the strategy and say which pairs will be hidden
  # and which ones will be reported to the mechanism.
  # 
  # Args:
  #   rke: An RKE object
  #   strategy.str : A strategy profile as string, "tttcc" etc..
  # Returns:
  #   A strategy object (see terminology)
  ret = list()
  if(strategy.str == "t") {
    ret$hide = c()
  } else if(strategy.str == "c") {
    # TODO(ptoulis): If the max matching is timed-out, this will return empty match
    # which is equivalent to being truthful. Needs a fix?
    matching = max.matching(rke, include.3way=include.3way)
    ret$hide = get.matching.ids(matching)
  }  else if(strategy.str == "r") {
    pairs.AB = rke.filter.pairs(rke, attr="desc", value="A-B")
    pairs.BA = rke.filter.pairs(rke, attr="desc", value="B-A")
    pairs.O = rke.filter.pairs(rke, attr="pair.type", value="O")
    ## want to hide the "short" side of R
    hide.R = pairs.AB
    if(length(pairs.BA) < length(pairs.AB))
      hide.R = pairs.BA
    ret$hide = c(hide.R)
  } 
  all.pairs = rke.pair.ids(rke)
  ret$report = setdiff(all.pairs, ret$hide)
  CHECK_strategy(ret, all.pairs)
  return(ret)
}

play.strategies = function(rke.list, strategy.str,
                           include.3way=F,
                           verbose=F) {
  # Initializes the RKE lists before the mechanism is run.
  # Reads the strategies and then plays them on each element of rke.list.
  # Calls play.strategy()
  #
  # Returns:
  #   A list(hospital.id => strategy)
  CHECK_rke.list(rke.list)
  strategies = strsplit(strategy.str,split="")[[1]]
  CHECK_MEMBER(strategies, c("t","c", "r", "b"), msg="Correct strategy spec")
  logthis("Strategies", verbose)
  logthis(strategies, verbose)
  hids <- get.rke.list.hospital.ids(rke.list)
  logthis("Hospital ids are: ", verbose)
  logthis(hids, verbose)
  strategy.list = llply(hids, function(h) play.strategy(rke.list[[h]],
                                                    strategy.str=strategies[h],
                                                    include.3way=include.3way))
  return(strategy.list)
}

get.hospitals.utility <- function(rke.all, matched.ids) {
  # Calculates utility per hospital given the supplied matched pair ids.
  # Returns a m x 1 vector of utilities
  matched.hospitals = rke.pairs.hospitals(rke.all, matched.ids)
  freq = sapply(rke.hospital.ids(rke.all), function(hid) length(which(matched.hospitals==hid)))
  return(as.matrix(freq, ncol=1))
}

##  Runs a mechanism
Run.Mechanism = function(kpd, mech, include.3way=F, verbose=F) {
  logthis(sprintf("Mechanism %s. Checking KPD", mech), verbose)
  CHECK_kpd(kpd)
  
  # unload
  reported.rke.list = kpd$reported.pool$rke.list
  reported.rke.all = kpd$reported.pool$rke.all
  rke.list = kpd$real.pool$rke.list
  rke.all = kpd$real.pool$rke.all
  
  ## 3. Run the mechanism
  logthis("Unload complete. Running the mechanism", verbose)
  matching = do.call(mech, args=list(rke.pool=kpd$reported.pool))
  logthis(sprintf("Matching status %s", matching$status), verbose)
  if (get.matching.status(matching) != "OK") {
    print(sprintf("Error happened in mechanism %s. Saving", mech))
    save(kpd, file="debug/kpd-%s.Rdata", kpd)
    stop("Quitting.")
  }
  matched.pairs = get.matching.ids(matching)
  ## Sometimes we get a list instead of a vector.
  logthis("Checking matched pairs", verbose)
  mech.out.ids = matched.pairs
  
  ## 4. Compute utility from mechanism
  logthis("Computing utilities", verbose)
  Util = get.hospitals.utility(reported.rke.all, mech.out.ids)
  logthis("Computed hospital Utility:", verbose)
  logthis(Util, verbose)
  hids <- rke.hospital.ids(rke.all)
  ## 5.  Utility from final internal matches.   
  for(h in hids) {
    rke.h = rke.list[[h]]
    # 1. For hospital h, find which pairs were matched
    hosp.matched.ids = intersect(rke.hospital.pairs(reported.rke.all, h),
                                 mech.out.ids)
    # 2. Remove from hidden part
    rke.remainder = rke.remove.pairs(rke.h, hosp.matched.ids)
    # 3. Perform maximum matching on hidden part
    logthis(sprintf("Remainder graph has %d pairs", rke.size(rke.remainder)), 
            verbose)
    matching = max.matching(rke.remainder, include.3way=include.3way)
    utility = get.matching.utility(matching)
    logthis(sprintf("Extra utility for hospital %d = %d", h, utility), verbose)
    # 4. Update utilities
    Util[h] = Util[h] + utility
  }
  logthis("Finaly Utility:", verbose)
  logthis(Util, verbose)
  return(Util)
}

rCM <- function(rke.pool, include.3way=F) {
  # Implementation of rCM. Returns an array of matched pair ids.
  CHECK_rke.pool(rke.pool)
  rke.list <- rke.pool$rke.list
  rke.all <- rke.pool$rke.all
  ## 1. Simply calculate a maximum-matching (this will shuffle the edges by default)
  m.all =  max.matching(rke.all, include.3way=include.3way)
  return(m.all)
}


compute.ir.constraints = function(rke.list, pair.types=c()) {
  # Given the list of RKE objects, compute the IR constraints for each hospital
  # Returns empty data-frame if not types, or empty rke.list
  #
  # Args:
  #   pair.types is an array of "pair types" which we compute constraints for.
  #
  # Returns:
  #   A data-frame that has (pc, hospital, internal.matches, pair.type) as fields, where
  #   pc = the pair code
  #   hospital = the hospital #id
  #   internal = #internal matches for this hospital and pc.
  #   pc  pair.type desc  hospital  internal.matches
  #   1     S       O-O        2           0
  #   2     U       A-O       2           1         ...
  stop("NEED TO use new definitions.")
  m = length(rke.list)
  CHECK_rke.list(rke.list)  # will check whether the hospital ids are correct.
  out = data.frame(pc=c(), pair.type=c(), hospital=c())
  if(length(pair.types)==0 || m == 0)
    return(out)
  
  for(hid in 1:length(rke.list)) {
    # xCM mechanism starts
    #  1. Compute max matching internally in S
    rke.h = rke.list[[hid]]
    CHECK_SETEQ(rke.hospital.ids(rke.h), c(hid), msg="Correct hospital id")
    for (type in c("S", "R")) {
      if(type %in% pair.types) {
        sub.pairs = rke.filter.pairs(rke.h, attr="pair.type", value=type)
        subrke = rke.keep.pairs(rke.h, pair.ids=sub.pairs)
        match = max.matching(subrke)
        CHECK_TRUE(match$results == "OK", msg="Matching should be performed")
        match = match$match  # ah well.
        tab.pc = as.matrix(table(match$pc))
        # data.frame of pc frequencies
        tab.pc = data.frame(pc=as.numeric(rownames(tab.pc)), freq=tab.pc[ , 1])
        CHECK_UNIQUE(tab.pc$pc)
        num.results = nrow(tab.pc)
        if (num.results > 0)
          out = rbind(out, list(pc=tab.pc$pc,
                                pair.type=rep(type, num.results),
                                desc=pc.to.desc(tab.pc$pc),
                                hospital=rep(hid, num.results),
                                internal.matches=tab.pc$freq))
      }
    }
  }## for every hospital
  return(out)
}

g.share = function(demand, supply) {
  # Fair sharing of supply.
  # The algorithm allocates supply fairly to each hospital.
  # Computes and allocation i.e. how many #items each hospital gets.
  # sum(allocation) = supply has to be.
  #
  # Args:
  #   z = m x 1 demand from agents
  #   x = supply (integer)
  #
  # Returns: 
  #   m x 1 vector of allocation
  # TO-DO(ptoulis): Unit tests for this one?
  #
  # J = { those ids that have positive demand }
  set.J = which(demand > 0)
  # allocation vector
  alloc = rep(0, length(demand))
  while(supply > 0 && length(set.J) > 0) {
    if(supply >= length(set.J)) {
      alloc[set.J] = alloc[set.J]+1
      supply = supply - length(set.J)
      set.J = setdiff(set.J, which(alloc >= demand))
    } else {
      j.subset = sample(set.J, size=supply, replace=F)
      alloc[j.subset] = alloc[j.subset]+1
      supply = 0
    }
  }
  return(alloc)
}

## TO-DO(ptoulis): Don't really care much about correctness
## as long as it has good properties, we can treat it as a black box.
xCM <- function(rke.pool) {
  # unload
  rke.list = rke.pool$rke.list
  rke.all = rke.pool$rke.all
  m = length(rke.list)
  matched.all.ids = c()
  ##  1. Compute IR constraints
  ir.constraints = compute.ir.constraints(rke.list, types=c("S", "R"))
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

  match.r = list()
  q = 0    
  ## Pair code (useful when setting constraints)
  pc.AB = pair.code(list(donor="A", patient="B"))
  pc.BA = pair.code(list(donor="B", patient="A"))
  pc.R = c(pc.AB, pc.BA)
  ## TO-DO(ptoulis): Bug empty graphs cycle forever
  keep.running = T
  while(keep.running) {
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
    keep.running = "INFEASIBLE" %in% names(match.r);
    q = q + 1
  }
  
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
    remove.edges = get.external.edges(rke.all, Gh.remainder)
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
  #warning("UD-lottery is not unit-tested.")
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
  #logwarn("Bonus is not unit-tested.")
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

