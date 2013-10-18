# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Contains code implementing different KPD mechanisms. 
# Currently supported:   rCM,  xCM
# A mechanism receives a combined donor-patient graph and outputs a matching. 
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
  #   A strategy object (see terminology) -- has "report" and "hide" fields
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
  hids <- rke.list.hospital.ids(rke.list)
  logthis("Hospital ids are: ", verbose)
  logthis(hids, verbose)
  strategy.list = llply(hids, function(h) play.strategy(rke.list[[h]],
                                                    strategy.str=strategies[h],
                                                    include.3way=include.3way))
  return(strategy.list)
}

get.hospitals.utility <- function(rke.all, matched.ids, all.hospital.ids) {
  # Calculates utility per hospital given the supplied matched pair ids.
  # Returns a m x 1 vector of utilities
  matched.hospitals = rke.pairs.hospitals(rke.all, matched.ids)
  freq = sapply(all.hospital.ids, function(hid) length(which(matched.hospitals==hid)))
  return(as.matrix(freq, ncol=1))
}

##  Runs a mechanism
Run.Mechanism = function(kpd, mech, include.3way, verbose=F) {
  # Runs a mechanism for a specified KPD
  # Returns:
  #     A matching object.
  warning("No unit tests for RunMechanism")
  logthis(sprintf("Mechanism %s. Checking KPD", mech), verbose)
  CHECK_kpd(kpd)
  
  # Init total matching.
  total.matching <- empty.match.result(empty.rke())
  
  # unload
  reported.rke.list = kpd$reported.pool$rke.list
  reported.rke.all = kpd$reported.pool$rke.all
  rke.list = kpd$real.pool$rke.list
  rke.all = kpd$real.pool$rke.all
  
  ## 3. Run the mechanism. Get the matching
  logthis("Unload complete. Running the mechanism", verbose)
  mech.matching = do.call(mech, args=list(rke.pool=kpd$reported.pool,
                                     include.3way=include.3way))
  
  total.matching = add.matching(total.matching, mech.matching)
  
  # Process matching.
  if (get.matching.status(mech.matching) != "OK") {
    print(sprintf("Error happened in mechanism %s. Saving", mech))
    save(kpd, file="debug/kpd-%s.Rdata", kpd)
    stop("Quitting.")
  }
  
  matched.pairs = get.matching.ids(mech.matching)
  ## Sometimes we get a list instead of a vector.
  logthis("Checking matched pairs", verbose)
  mech.out.ids = matched.pairs
  
  ## 4. Compute utility from mechanism
  logthis("Computing utilities", verbose)
  all.hospitals <- rke.list.hospital.ids(rke.list)
  hids <- all.hospitals # important to have all hospitals
  ## 5.  Utility from final internal matches.   
  for(h in hids) {
    rke.h = rke.list[[h]]
    # 1. For hospital h, find which pairs were matched
    hosp.matched.ids = intersect(rke.hospital.pair.ids(reported.rke.all, h),
                                 mech.out.ids)
    # 2. Remove from hidden part
    rke.remainder = rke.remove.pairs(rke.h, hosp.matched.ids)
    # 3. Perform maximum matching on hidden part
    logthis(sprintf("Remainder graph has %d pairs", rke.size(rke.remainder)), 
            verbose)
    matching = max.matching(rke.remainder, include.3way=include.3way)
    
    # 0.2 Add this matching to the total
    total.matching <- add.matching(total.matching, matching)
    
  }
  return(total.matching)
}

rCM <- function(rke.pool, include.3way=F) {
  # Implementation of rCM. Returns an array of matched pair ids.
  # loginfo(sprintf("Running xCM 3-chain=%s", include.3way))
  CHECK_rke.pool(rke.pool)
  rke.all <- rke.pool$rke.all
  
  # Total matching computed by the mechanism.
  total.matching = empty.match.result(empty.rke())
  
  ## 1. Simply calculate a maximum-matching (this will shuffle the edges by default)
  m.all =  max.matching(rke.all, include.3way=include.3way)
  
  ## 0.1 Add matching
  total.matching <- add.matching(total.matching, m.all)
  
  return(total.matching)
}


compute.ir.constraints = function(rke.pool, pair.types=c(), include.3way=F) {
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
  #   2     U       A-O        2           1         ...
  CHECK_rke.pool(rke.pool)  # will check whether the hospital ids are correct.
  if (length(pair.types) == 0)
    warning("No pair types defined. Returning empty constraints.")
  rke.list = rke.pool$rke.list
  m = length(rke.list)
  out = data.frame(pc=c(), pair.type=c(), hospital=c(), internal.matches=c())
  if(length(pair.types)==0 || m==0) {
    warning("Empty RKE list.")
    return(out)
  }
  for(hid in 1:length(rke.list)) {
    # 1. Compute max matching internally in S
    rke.h = rke.list[[hid]]
    for (type in intersect(pair.types, c("S", "R"))) {
      # Take the pairs the belong to this group
      sub.pair.ids = rke.filter.pairs(rke.h, attr="pair.type", value=type)
      subrke = rke.keep.pairs(rke.h, pair.ids=sub.pair.ids)
      matching = max.matching(subrke, include.3way=include.3way)
      CHECK_TRUE(get.matching.status(matching) == "OK", msg="Matching OK?")
      match = matching$match
      tab.pc = as.matrix(table(match$pc))  # tabulate the PC codes.
      # data.frame of pc frequencies
      tab.pc = data.frame(pc=as.numeric(rownames(tab.pc)), freq=tab.pc[ , 1])
      CHECK_UNIQUE(tab.pc$pc)
      num.results = nrow(tab.pc)
      if (num.results > 0)
        out = rbind(out, list(pc=tab.pc$pc,
                              hospital=rep(hid, num.results),
                              internal.matches=tab.pc$freq))
    }
  }## for every hospital
  rownames(out) <- NULL
  out$desc = sapply(out$pc, pc.to.desc)
  out$pair.type = pc.to.pair.type(out$pc)
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
  #
  # J = { those ids that have positive demand }
  warning("Not unit-tests for g share")
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

compute.Rsubgraph.constraints <- function(ir.constraints, rke.pool) {
  # Constraints for the R subgraph
  # Returns:
  #   pc  hospital  internal.matches  pair.type  y
  #
  # pc = pair code, hospital = hospital id
  # internal.matches = #matched in R subgraph
  # Should be pair.type == R
  # y = Allocation after running the lottery for the unmatched R pairs.
  m = length(rke.pool$rke.list)
  empty.cons <- function() {
    x <- expand.grid(hospital=c(1:m), pc=c(7,10))
    x$pair.type <- "R"
    x$desc <- pc.to.desc(pcs=x$pc)
    x$internal.matches <- 0
    x$y <- 0
    return(x)
  }
  cons <- empty.cons()
  if(nrow(ir.constraints) == 0) {
    cons <- empty.cons()
  } else {
    cons <- subset(ir.constraints, pair.type=="R")
  }
  tmp.frame <- expand.grid(hospital=c(1:m), pc=c(7,10))
  cons = join(cons, tmp.frame, type="right", by=c("pc", "hospital"))
  cons$internal.matches[is.na(cons$internal.matches)] <- 0
  cons$pair.type <- "R"
  cons$desc <- pc.to.desc(pcs=cons$pc)
  count.h.pairs <- function(hid, pc) {
    length(rke.filter.pairs(rke.pool$rke.all, attrs=c("hospital", "pc"), values=c(hid, pc)))
  }
  CHECK_EQ(nrow(cons), 2 * m, msg="Correct no. of AB/BA #constraints")
  cons$unmatched <- sapply(1:nrow(cons), function(i) with(cons[i, ], count.h.pairs(hospital, pc)))
  cons$unmatched <- cons$unmatched - cons$internal.matches
  CHECK_GE(cons$unmatched, 0, msg="Non-negative #matches")
  consAB <- subset(cons, desc=="A-B")
  consBA <- subset(cons, desc=="B-A")
  
  xAB <- sum(consAB$unmatched)
  xBA <- sum(consBA$unmatched)
  yAB <- rep(0, m)
  yBA <- rep(0, m)
  if (xAB >= xBA)
    yAB <- g.share(demand=consAB$unmatched, supply=xBA)
  else
    yBA <- g.share(demand=consBA$unmatched, supply=xAB)
  consAB$y <- yAB
  consBA$y <- yBA
  cons <- rbind(consAB, consBA) 
  
  return(cons)
}

## TO-DO(ptoulis): Don't really care much about correctness
## as long as it has good properties, we can treat it as a black box.
#
# Returns: LIST(matching=total.matching)
xCM <- function(rke.pool, include.3way=F, verbose=F) {
  CHECK_rke.pool(rke.pool)
  # loginfo(sprintf("Running xCM 3-chain=%s", include.3way))
  # unload
  rke.list = rke.pool$rke.list
  rke.all = rke.pool$rke.all
  
  # Total matching computed by the mechanism.
  total.matching = empty.match.result(empty.rke())
  
  m = length(rke.list)  # no. of hospitals
  ##  1. Compute IR constraints
  ir.constraints = compute.ir.constraints(rke.pool, pair.types=c("S", "R"), include.3way=include.3way)
  
  # 2.  Match S internally
  s.subrke = rke.subgraph(rke.all, pair.type="S")
  s.constraints = subset(ir.constraints, pair.type=="S")
  s.matching = max.matching(s.subrke,
                            ir.constraints=s.constraints,
                            include.3way=include.3way)
  # 0.1 Add matching information from S
  total.matching <- add.matching(total.matching, s.matching)
  
  ## 3.   Match R internally
  r.subrke = rke.subgraph(rke.all, pair.type="R")
  r.constraints = compute.Rsubgraph.constraints(ir.constraints, rke.pool=rke.pool)
  # R-constraints have additionally the "y" lottery allocation (see notes)
  CHECK_EQ(nrow(r.constraints), 2 * m, msg="AB and BA pairs for each hospital")
  r.matching = list()
  q = 0    
  loop.ended = F
  while(!loop.ended) {
    loop.r.constraints = r.constraints
    loop.r.constraints$unmatched <- r.constraints$unmatched + max(0, r.constraints$y - q)
    r.matching = max.matching(r.subrke, ir.constraints = loop.r.constraints,
                              include.3way=include.3way)
    loop.ended = get.matching.status(r.matching) == "OK"
    q = q + 1
  }
  # 0.2 Add matching information from R
  total.matching <- add.matching(total.matching, r.matching)
  ## remove some stuff that are not needed anymore
  rm(ir.constraints)
  logthis(sprintf("Total S matches = %d", length(get.matching.ids(s.matching))), verbose)
  logthis(sprintf("Total R matches = %d", length(get.matching.ids(r.matching))), verbose)

  #  4. Almost regular matching to the remainder
  ## Match OD's individually.
  for(hid in 1:m) {
    matched.hospital.ids <- intersect(total.matching$match$pair.id,
                                      rke.hospital.pair.ids(rke.all, hid))
    rke.h <- rke.remove.pairs(rke=rke.list[[hid]],
                              rm.pair.ids=matched.hospital.ids)
    internal.matching = max.matching(rke.h, 
                                     include.3way=include.3way,
                                     regular.matching=T)
    # 0.3 Add matching information from internal matching.
    total.matching <- add.matching(total.matching, internal.matching)
  }
  
  return(matching=total.matching)
}


ud.lottery = function(rke, ud.pc, 
                      Hn, theta, QS) {
  # Under-demanded lottery as in Ashlagi-Roth (2013)
  #
  # Args
  #  ud.pc = a Pair code of U pair
  #  Hn = array of hospital ids
  #  theta = supply
  #  QS = the QS data structure (see Bonus.QS). This holds information 
  #     of U pairs that hospitals have or can match
  #
  # Returns
  #   An array of U pair ids (that belong to hospitals in Hn)
  #
  # warning("UD-lottery is not unit-tested.")

  # X-Y pairs in rke for the Hn hospitals.
  rke.XY.pair.ids = subset(rke$pairs, pc==ud.pc & hospital %in% Hn)$pair.id
  # Bh.XY = LIST (per hospital).
  # Bh.XY[[hid]] = { set of XY pairs for hospital hid }
  Bh.XY = list()
  # Sh.XY = LIST (per hospital).
  # Sh.XY[[hid]] = { set of XY pairs that are matched internally }
  Sh.XY = list()
  ## Contains # balls /hospital     [1,1,1,3,3,3,3,3,5,5,5, ...]   NOTE: this is index by a hospital in Hn
  set.J = c()
  ## Initialize the Sh.XY, Bh.XY
  for(hid in Hn) {
    Sh.XY[[hid]] = QS$S[[hid]][[ud.pc]]  # set of matched XY pairs.
    set.J = c(set.J, rep(hid, QS$Q[[hid]][ud.pc]))  # add balls for hospital
    Bh.XY[[hid]] = rke.filter.pairs(rke, attrs=c("pc", "hospital"),
                                    values=c(ud.pc, hid))
    ## Sanity check (matched ids subset of all ids)
    CHECK_MEMBER(Sh.XY[[hid]], Bh.XY[[hid]])
    CHECK_EQ(length(Bh.XY[[hid]]), QS$Q[[hid]][ud.pc])
  }
  # Sum of all candidate matched XY pairs across hospitals
  get.sum = function() sum(laply(Sh.XY, length))
  plateau.count = 0
  # Main loop of lottery
  # Do not remove balls from J (that is the definition)
  # J = 1,1,1,2,2,2,3,3,4,...  (contains hospital ids)
  if (length(set.J) > 0)
    while(get.sum() < theta && plateau.count < 10) {
      h.sample = sample(set.J, size=1)
      # set.J = set.J[-which(set.J==h.sample)[1]]
      CHECK_MEMBER(Sh.XY[[h.sample]], Bh.XY[[h.sample]],
                   msg="Matched ids subset of all")
      
      avail = setdiff(Bh.XY[[h.sample]], Sh.XY[[h.sample]])
      if(length(avail) > 0) {
        xy.pair = uniform.sample(avail)
        Sh.XY[[h.sample]] = c(Sh.XY[[h.sample]], xy.pair)
        plateau.count=0
      } else {
        plateau.count = plateau.count +1 
      }
    }
  
  ret.pair.ids = c()
  # Return the union.
  for(hid in Hn) {
    ret.pair.ids = c(ret.pair.ids, Sh.XY[[hid]])
  }
  CHECK_MEMBER(ret.pair.ids, rke.XY.pair.ids, msg="Candidate matched subset of all")
  return(ret.pair.ids)
}

Bonus.QS = function(rke.pool, include.3way=F)  {
  # Compute Qh and Sh sets internally
  # 
  # Returns a LIST that has:
  #   Q[hid][X-Y] = how many #X-Y  in hospital hid 
  #   S[hid][X-Y] = { ids } of type X-Y which are included in a max regular matching.
  Qh = list()
  Sh = list()
  CHECK_rke.pool(rke.pool)
  # unload
  rke.list = rke.pool$rke.list
  rke.all = rke.pool$rke.all
  all.hospital.ids = rke.list.hospital.ids(rke.list)
  CHECK_EQ(all.hospital.ids, 1:max(unique(all.hospital.ids)))
  # pair codes of Underdemanded pairs (see terminology)
  ud.pcs = subset(kPairs, pair.type=="U")$pc
  
  for(hid in all.hospital.ids) {
    Qh[[hid]] = rep(0, length(kPairCodes))
    Sh[[hid]] = list()
    hospital.pair.ids = rke.hospital.pair.ids(rke.all, hospital.id=hid)
    hospital.pcs = subset(rke.all$pairs, hospital==hid)$pc
    # 2 . Compute a regular matching on the specific hospital (internal)
    internal.matching = max.matching(rke.list[[hid]], include.3way=include.3way,
                                     regular.matching=T)
    ## Iterate over all UD codes
    for(pc in ud.pcs) {
      hpairs.of.type = hospital.pair.ids[which(hospital.pcs==pc)]
      Qh[[hid]][pc] = length(hpairs.of.type)
      Sh[[hid]][[pc]] = intersect(get.matching.ids(internal.matching),
                                  hpairs.of.type)
    }
  }
  ##  TO-DO(ptoulis): How do you know if you computed Qh, Sh correctly???
  return(list(Q=Qh, S=Sh))  
}

## Bonus mechanism. Ashlagi & Roth (2013)
Bonus = function(rke.pool, include.3way=F) {
  CHECK_rke.pool(rke.pool)
  # unload
  rke.list = rke.pool$rke.list
  rke.all = rke.pool$rke.all
  
  # Total matching computed by the mechanism.
  total.matching = empty.match.result(empty.rke())
    
  ## 0. Initialize mechanism
  m = length(rke.list)  # no. of hospitals.

  ## 1. Match S pairs internally
  s.subrke = rke.subgraph(rke.all, pair.type="S")
  s.matching = max.matching(s.subrke, include.3way=include.3way)
  
  ## 0.1 Add S-matching
  total.matching <- add.matching(total.matching, s.matching)
  
  ## 2. Match R pairs
  # Compute constraints.
  R.constraints = compute.ir.constraints(rke.pool, pair.types=c("R"))
  r.subrke = rke.subgraph(rke.all, pair.type="R")
  r.matching = max.matching(r.subrke, include.3way=include.3way, 
                            ir.constraints=R.constraints)
  ## 0.1 Add R-matching
  total.matching <- add.matching(total.matching, r.matching)
  
  #  3. Match OD/UD pairs.
  #   3a. Split hospital in two sets
  k = as.integer(m/2)
  H.sets = list()
  H.sets[[1]] = sample(1:m, size=k, replace=F)
  H.sets[[2]] = setdiff(1:m, H.sets[[1]])

  # All Underdemanded PC codes
  ud.pcs = subset(kPairs, pair.type=="U")$pc
  # Compute Q-S:
  #   Q[hid][X-Y] = how many #X-Y in hospital hid
  #   S[hid][X-Y] = { ids of matched X-Y for hid }
  QS.obj = Bonus.QS(rke.pool=rke.pool, include.3way=include.3way)
  ## For all under-demanded pairs  X-Y
  for(i in ud.pcs) {
    # i is a PC (pair code)
    # Note: follow the notation in Ashlagi&Roth about Bonus
    # Y-X (O pairs) and X-Y (U pairs) are reciprocal.
    all.XY = rke.filter.pairs(rke.all, attrs=c("pc"), values=c(i))
    all.YX = rke.filter.pairs(rke.all, attrs=c("pc"), values=pc.reciprocal(i))
    # Iterate through hospital groups
    for(j in 1:2) {
      # Get YX pairs that belong to the H2
      # and XY pairs that belong to H1 (indexed by "j")
      # tau.Bhj.XY = { ids that belong to H1 and are of type XY }
      tau.BHj.XY = subset(rke.all$pairs,
                          hospital %in% H.sets[[j]] & pair.id %in% all.XY)$pair.id
      # tau.BHother.XY = { ids that belong to H2 and are of type YX }
      tau.BHother.YX = subset(rke.all$pairs,
                              hospital %in% H.sets[[3-j]] & pair.id %in% all.YX)$pair.id
      
       # no. of such pairs (supply of YX pairs to be used in lottery)
      theta.j.YX = length(tau.BHother.YX)
      
      ## Run the under-demanded lottery  (contains pair ids)
      S.XY = ud.lottery(rke=rke.all,
                        ud.pc=i,
                        Hn=H.sets[[j]],
                        theta=theta.j.YX,
                        QS=QS.obj)
      
      CHECK_DISJOINT(S.XY, tau.BHother.YX, msg="Not even the same type")
      CHECK_MEMBER(S.XY, tau.BHj.XY, msg="Matched subset of all")
      ## the only ids to be considered in the matching
      XYYX.ids = union(S.XY, tau.BHother.YX)
      # Here you match   X-Y pairs from Hj   with   Y-X pairs from Hother
      subrke = rke.keep.pairs(rke.all, pair.ids=XYYX.ids)
      Mj = max.matching(subrke, include.3way=include.3way, regular.matching=F)
      
      ## 0.3 Add new matching
      total.matching <- add.matching(total.matching, Mj)
    }
  }
  return(total.matching)
}