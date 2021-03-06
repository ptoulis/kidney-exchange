# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
#
# Contains code implementing different KPD mechanisms. 
# Currently supported:   rCM,  xCM
# A mechanism receives a combined donor-patient graph and outputs a matching. 

# Make sure these are loaded.
# source("rke.R")
# source("matching.R")

kImplementedKPDMechanisms = c("rCM", "selfCM", "xCM", "Bonus")

kpd.create <- function(rke.pool, strategy.str, include.3way, verbose=F) {
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
  hospital.ids <- rke.list.hospital.ids(rke.list)
  CHECK_EQ(nchar(strategy.str), m, msg="Correct #of strategies")
  CHECK_SETEQ(1:m, hospital.ids)
  logfine("Playing strategies")
  
  hospital.action = play.strategies(rke.list=rke.list, 
                                    strategy.str=strategy.str, 
                                    include.3way=include.3way)
  
  reported.rke.list = list()
  ## hidden.pairs = ids in terms of rke.all  that were hidden
  all.hidden.pairs = c()
  #  1b.  Populate the reported and hidden graphs.
  for(hid in hospital.ids) {
    hospitalAction = hospital.action[[hid]]
    logfine(sprintf("Hospital %d is hiding %d pairs", hid, length(hospitalAction$hide)))
    reported.rke.list[[hid]] = rke.remove.pairs(rke.list[[hid]], hospitalAction$hide)
    CHECK_DISJOINT(all.hidden.pairs, hospitalAction$hide)
    all.hidden.pairs = c(all.hidden.pairs, hospitalAction$hide)
  }
  ##  Save the "kpd" object
  ## This is important in order to compare mechanisms
  ## on the same sets of hospitals and pooled graphs.
  kpd = list()
  kpd$reported.pool = list(rke.list=reported.rke.list,
                           rke.all=rke.remove.pairs(rke.all, all.hidden.pairs))
  kpd$real.pool = list(rke.list=rke.list, rke.all=rke.all)
  CHECK_kpd(kpd)
  return(kpd)
}

play.strategy <- function(rke, strategy.str, include.3way) {
  # For the RKE, play the strategy and say which pairs will be hidden
  # and which ones will be reported to the mechanism.
  # 
  # Args:
  #   rke: An RKE object
  #   strategy.str : A strategy profile as string, "tttcc" etc..
  # Returns:
  #   A strategy object =LIST(hide, report) (see terminology) -- has "report" and "hide" fields
  #      hide, report = ARRAY(pair ids)
  CHECK_MEMBER(strategy.str, c("t", "c", "r"))
  ret = list()
  if(strategy.str == "t") {
    ret$hide = c()
  } else if(strategy.str == "c") {
    # TODO(ptoulis): If the max matching is timed-out, this will return empty match
    # which is equivalent to being truthful. Needs a fix?
    # Update (12/2013): Canonical should be also be regular.
    matching = max.matching(rke, include.3way=include.3way, regular.matching=T)
    ret$hide = get.matching.ids(matching)
  }  else if(strategy.str == "r") {
    pairs.AB = rke.filter.pairs(rke, attr="desc", value="A-B")
    pairs.BA = rke.filter.pairs(rke, attr="desc", value="B-A")
    # want to hide the "short" side of R
    hide.R = pairs.AB
    if(length(pairs.BA) < length(pairs.AB))
      hide.R = pairs.BA
    ret$hide = c(hide.R)
  } 
  all.pair.ids = rke.pair.ids(rke)
  ret$report = setdiff(all.pair.ids, ret$hide)
  CHECK_strategy(ret, all.pair.ids)
  return(ret)
}

play.strategies = function(rke.list, strategy.str, include.3way, verbose=F) {
  # Initializes the RKE lists before the mechanism is run.
  # Reads the strategies and then plays them on each element of rke.list.
  # Calls play.strategy()
  #
  # Returns:
  #   A LIST(hospital.id => strategy)
  CHECK_rke.list(rke.list)
  strategies = strsplit(strategy.str,split="")[[1]]
  CHECK_MEMBER(strategies, c("t", "c", "r"), msg="Correct strategy spec")
  CHECK_EQ(length(strategies), length(rke.list), msg="#strategies = #hospitals")
  hids <- rke.list.hospital.ids(rke.list)
  strategy.list = llply(hids, function(h) {
    play.strategy(rke.list[[h]],
                  strategy.str=strategies[h],
                  include.3way=include.3way)
  })
  return(strategy.list)
}

Run.Mechanism = function(kpd, mech, include.3way, verbose=F) {
  # Runs a mechanism for a specified KPD
  # Returns:
  #   A LIST of MATCHING objects:
  #     LIST(total.matching, mech.matching, internal.matchings)
  #   total.matching = overall MATCHING from running the mechanism
  #   mech.matching = MATCHING of the mechanism
  #   internal.matchings = LIST of m MATCHING elements, one for each hospital
  #     it is the internal matchings that hospitals make on the remainders.
  logfine("Checking KPD object")
  CHECK_kpd(kpd)
  
  # Init total matching.
  total.matching <- empty.match.result(empty.rke())
  mech.matching <- empty.match.result(empty.rke())
  internal.matchings <- list()
  
  # unload
  reported.rke.list = kpd$reported.pool$rke.list
  reported.rke.all = kpd$reported.pool$rke.all
  rke.list = kpd$real.pool$rke.list
  rke.all = kpd$real.pool$rke.all
  logfine("Compute mechanism matching.")
  ## 3. Run the mechanism. Get the matching
  mech.matching = do.call(mech, args=list(pool=kpd$reported.pool,
                                     include.3way=include.3way))
  
  total.matching = add.matching(total.matching, mech.matching)
  logfine(sprintf("Mech matching computed. Total %d matches", nrow(total.matching$match)))
  
  # Zero-tolerance to errors.
  if (get.matching.status(mech.matching) != "OK") {
    print(sprintf("Error happened in mechanism %s. Saving", mech))
    save(kpd, file="debug/kpd-%s.Rdata", kpd)
    stop("Quitting.")
  }
  
  ## 4. Compute utility from mechanism
  hids <- rke.list.hospital.ids(rke.list)
  matched.pair.ids = get.matching.ids(total.matching)
  
  ## 5.  Utility from final internal matches. Recourse step
  for(h in hids) {
    rke.h = rke.list[[h]]
    # 1. For hospital h, find which pairs were matched
    hosp.matched.ids = intersect(rke.pair.ids(rke.h), matched.pair.ids)
    # 2. Remove from hidden part
    rke.remainder = rke.remove.pairs(rke.h, rm.pair.ids=hosp.matched.ids)
    logfine(sprintf("RunMechanism-Recourse: Matching for hospital %d, total %d pairs remaining", h, rke.size(rke.remainder)))
    # 3. Perform maximum matching on hidden part
    matching = max.matching(rke.remainder,
                            include.3way=include.3way,
                            regular.matching=T)
    internal.matchings[[h]] <- matching
    # 0.2 Add this matching to the total
    total.matching <- add.matching(total.matching, matching)
  }
  
  if(verbose) {
    print("Mech matching information")
    print(mech.matching$information)
    print("Total matching information")
    print(total.matching$information)
    print("Breakdown by hospital")
    for(h in hids) {
      print(sprintf("Hospital %d", h))
      print(internal.matchings[[h]]$information)
    }
  }
  
  return(list(total.matching=total.matching,
              mech.matching=mech.matching,
              internal.matchings=internal.matchings))
}

rCM <- function(pool, include.3way=F) {
  # Implementation of rCM. 
  # Returns a MATCHING object. (all mechanisms should return a matching object)
  CHECK_rke.pool(pool)
  rke.all <- pool$rke.all
  
  # total.matching = MATCHING object to be returned.
  total.matching = empty.match.result(empty.rke())
  ## 1. Simply calculate a maximum-matching (this will shuffle the edges by default)
  m =  max.matching(rke.all, include.3way=include.3way, regular.matching=T)
  
  total.matching = m  # a bit bogus, but good to have structure.
  return(total.matching)
}

selfCM <- function(pool, include.3way=F) {
  # Implementation of selfish-Centralized mechanism.
  #
  # Matches all hospitals internally (no pooling)
  CHECK_rke.pool(pool)
  rke.all <- pool$rke.all
  rke.list <- pool$rke.list
  # total.matching = MATCHING object to be returned.
  total.matching = empty.match.result(empty.rke())
  
  ## 1. Simply calculate a maximum-matching for each hospital internally.
  for(h in rke.list.hospital.ids(rke.list)) {
    rke.h = rke.list[[h]]
    m = max.matching(rke.h, include.3way=include.3way, regular.matching=T)
    total.matching <- add.matching(total.matching, m)
  }
  
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
  
  hospital.ids <- rke.list.hospital.ids(rke.list)
  all.types = as.character(unique(kPairs$pair.type))
  CHECK_MEMBER(pair.types, all.types)
  for(hid in hospital.ids) {
    # 1. Compute max matching internally in S
    rke.h = rke.list[[hid]]
    for (type in pair.types) {
      # Take the pairs the belong to this group
      subrke = rke.subgraph(rke.h, pair.type=type)
      matching = max.matching(subrke, include.3way=include.3way)
      CHECK_TRUE(get.matching.status(matching) == "OK", msg="Matching OK?")
      # get the pair codes for this specific pair type.
      all.type.pcs = subset(kPairs, pair.type==type)$pc
      for(type.pc in all.type.pcs) {
        out = rbind(out, list(pc=type.pc,
                              hospital=hid,
                              internal.matches=length(which(matching$match$pc==type.pc))))
      }
    }  ## constraints for every hospital
  }
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
  #   pc  hospital  internal.matches  pair.type  unmatched y
  #
  # pc = pair code, hospital = hospital id
  # internal.matches = #matched in R subgraph
  # Should be pair.type == R
  # y = Allocation after running the lottery for the unmatched R pairs.
  m = length(rke.pool$rke.list)
  all.hids = rke.list.hospital.ids(rke.pool$rke.list)
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
  tmp.frame <- expand.grid(hospital=all.hids,
                           pc=subset(kPairs, pair.type=="R")$pc)
  cons = join(cons, tmp.frame, type="right", by=c("pc", "hospital"))
  cons$internal.matches[is.na(cons$internal.matches)] <- 0
  cons$pair.type <- "R"
  cons$desc <- pc.to.desc(pcs=cons$pc)
  CHECK_EQ(nrow(cons), 2 * m, msg="Correct no. of AB/BA #constraints")
  # count total R pairs for each hospital
  # Vector of size (m x 1)
  totalR <- sapply(1:nrow(cons), 
                   function(i) {
                     hid = cons[i, ]$hospital
                     hid.pc = cons[i, ]$pc
                     nrow(subset(rke.pool$rke.all$pairs, hospital==hid & pc==hid.pc))
                   })

  cons$unmatched <- totalR - cons$internal.matches
  CHECK_GE(cons$unmatched, 0, msg="Non-negative #matches")
  # consAB = zAB and consBA = zBA (but data-frames!)
  consAB <- subset(cons, desc=="A-B")
  consBA <- subset(cons, desc=="B-A")
  zAB = consAB$unmatched
  zBA = consBA$unmatched
  
  # New tweaks (30/12)
  # nAB = counts of AB pairs (over hospitals)
  nAB = sapply(all.hids, function(hid) {
    rke.h = rke.pool$rke.list[[hid]]
    nrow(subset(rke.h$pairs, desc=="A-B"))
  })
  nBA = sapply(all.hids, 
               function(hid) {
                 rke.h = rke.pool$rke.list[[hid]]
                 nrow(subset(rke.h$pairs, desc=="B-A"))
               })
  CHECK_TRUE(all(nAB + nBA >= cons$unmatched))
  # Ideal-matches: Assume perfect-matching short to long side.
  mh.ideal = sapply(1:length(all.hids), function(h) min(nAB[h], nBA[h]))
  # Real-matches: Just take the info from the constraints.
  m.real = sapply(all.hids, 
                  function(hid) {
                    tmp = subset(cons, hospital==hid & desc=="A-B")
                    if(nrow(tmp)==0) return(0)
                    return(tmp$internal.matches)
                  })
  zAB = nAB - mh.ideal
  zBA = nBA - mh.ideal
  delta = mh.ideal - m.real
  CHECK_TRUE(all(delta >= 0))
  CHECK_TRUE(all(zAB >= 0))
  CHECK_TRUE(all(zBA >= 0))
  ##  New tweaks end here.
  
  xAB <- sum(zAB)
  xBA <- sum(zBA)
  yAB <- rep(0, m)
  yBA <- rep(0, m)
  if (sum(nAB) >= sum(nBA))
    yAB <- g.share(demand=zAB, supply=xBA)
  else
    yBA <- g.share(demand=zBA, supply=xAB)
  consAB$y <- yAB + delta
  consBA$y <- yBA + delta
  cons <- rbind(consAB, consBA) 
  
  return(cons)
}

# Implements xCM mechanism (Toulis & Parkes, 2013)
# Returns: LIST(matching=total.matching)
xCM <- function(pool, include.3way=F, verbose=F) {
  logfine("Checking rke.pool")
  CHECK_rke.pool(pool)
  if(include.3way) {
    return(xCM3(pool, verbose=verbose))
  }
  # unload
  rke.list = pool$rke.list
  rke.all = pool$rke.all
  
  # Total matching computed by the mechanism.
  total.matching = empty.match.result(empty.rke())
  
  m = length(rke.list)  # no. of hospitals
  ##  1. Compute IR constraints
  logfine("Compute IR constraints for S-, R-subgraphs")
  ir.constraints = compute.ir.constraints(pool, pair.types=c("S", "R"), include.3way=F)
  
  # 2.  Match S internally
  s.subrke = rke.subgraph(rke.all, pair.type="S")
  logfine(sprintf("Compute matching on S-subgraph. Total %d pairs, %d edges",
                  rke.size(s.subrke), sum(s.subrke$edges$can.donate)))
  s.constraints = subset(ir.constraints, pair.type=="S")
  s.matching = max.matching(s.subrke,
                            ir.constraints=s.constraints,
                            include.3way=F)
  
  # 0.1 Add matching information from S
  total.matching <- add.matching(total.matching, s.matching)
  logfine("xCM: S-matching computed.")
  ## 3.   Match R internally
  r.subrke = rke.subgraph(rke.all, pair.type="R")
  logfine(sprintf("Compute matching on R-subgraph with constraints. Total %d pairs.",
                  rke.size(r.subrke)))
  r.constraints = compute.Rsubgraph.constraints(ir.constraints, rke.pool=pool)
  # R-constraints have additionally the "y" lottery allocation (see notes)
  CHECK_EQ(nrow(r.constraints), 2 * m, msg="AB and BA pairs for each hospital")
  
  r.matching = list()
  q = 0    
  loop.ended = F
  while(!loop.ended) {
    loop.r.constraints = r.constraints
    # loop.r.constraints$unmatched <- r.constraints$unmatched + max(0, r.constraints$y - q)
    loop.r.constraints$internal.matches = 
      sapply(1:nrow(r.constraints), function(i) {
        r.constraints$internal.matches[i] + max(0, r.constraints$y[i] - q)
      })

    r.matching = max.matching(r.subrke, ir.constraints = loop.r.constraints, include.3way=F)
    loop.ended = get.matching.status(r.matching) == "OK"
    q = q + 1
    logfine(sprintf("xCM: In R-subgraph loop. q=%d", q))
  }

  # 0.2 Add matching information from R
  total.matching <- add.matching(total.matching, r.matching)
  ## remove some stuff that are not needed anymore
  rm(ir.constraints)
  logfine("xCM: R-matching computed.")
  
  #  4. Almost regular matching to the remainder
  ## Match OD's individually.
  hospital.ids <- rke.list.hospital.ids(rke.list)
  for(hid in hospital.ids) {
    matched.hospital.ids <- intersect(get.matching.ids(total.matching),
                                      rke.pair.ids(rke.list[[hid]]))
    rke.h <- rke.remove.pairs(rke=rke.list[[hid]],
                              rm.pair.ids=matched.hospital.ids)
    internal.matching = max.matching(rke.h, include.3way=F, regular.matching=T)
    # 0.3 Add matching information from internal matching.
    total.matching <- add.matching(total.matching, internal.matching)
    logfine(sprintf("xCM: Computed internal matching for hospital %d", hid))
  }
  
  # match remainder.
  rke.remainder = rke.remove.pairs(rke.all, rm.pair.ids=total.matching$match$pair.id)
  remainder.m = max.matching(rke.remainder, include.3way=F, regular.matching=T)
  total.matching <- add.matching(total.matching, remainder.m)
  
  return(total.matching)
}

OAB.matched.in.OUU <- function(rke) {
  # Computes the matches of O-AB to U-U pairs.
  # In new notation i.e. (patient, donor) this would be AB-O pairs to U, U
  #
  # Returns LIST(OAB=(pair ids of OAB pairs that can be matched), ,
  #              U = pair ids that can be matched in OUU
  #              matching = the MATCHING object (see terminology.R)
  # we are using the old notation
  OABpair.ids = subset(rke$pairs, desc=="O-AB")$pair.id
  Upair.ids = subset(rke$pairs, pair.type=="U")$pair.id
  # get the OU RKE object
  OU.rke = rke.keep.pairs(rke, pair.ids=c(OABpair.ids, Upair.ids))
  empty.result = list(matching=empty.match.result(empty.rke()),
                      OAB=c(), U=c())
  if(length(OABpair.ids) == 0) return(empty.result)
  #
  # Compute a max-matching (w/ 3cycles) where U pairs take 2x weight.
  # This will get as many OUU exchanges as possible.
  m = max.matching(OU.rke, include.3way=T)
  if(m$utility==0) {
    return(list(matching=empty.match.result(empty.rke()),
                OAB=c(), U=c()))
  }
  # take the matched cycles and remove the "type" field.
  m3 = subset(m$matched.cycles, select=-type)
  # iterate all 3-way exchanges, count #OAB pairs matched.
  OAB.matched = as.numeric(apply(m3, 1, function(x) length(intersect(OABpair.ids, x))))
  # iterate all exchanes
  totalU.matched = as.numeric(apply(m3, 1, function(x) length(intersect(Upair.ids, x))))
  ind = intersect(which(OAB.matched==1), which(totalU.matched==2))
  # ind = index of 3-way cycles in m3, that match O-AB in OUU
  # need to take the pair ids.
  all.OAB.matched = c()
  all.U.matched = c()
  for(i in ind) {
    x = as.numeric(m3[i, ])
    all.OAB.matched <- c(all.OAB.matched, intersect(OABpair.ids, x))
    all.U.matched <- c(all.U.matched, intersect(Upair.ids, x))
  }
  return(list(matching=matching.keep.pair.ids(OU.rke, m, c(all.OAB.matched, all.U.matched)),
              OAB=all.OAB.matched,
              U=all.U.matched))         
}

xCM3 <- function(pool, verbose=F) {
  # Runs xCM with 3-way exchanges.
  #
  # The most important object is "total.matching" for every mechanism
  # In each step this is updated with all matches performed so far.
  # This is a MATCHING object (se terminology)
  CHECK_rke.pool(pool)
  # unload the RKE pool.
  rke.list = pool$rke.list  # LIST of individual RKE objects
  rke.all = pool$rke.all  # global RKE (pooled graph)
  
  # Total matching computed by the mechanism. This will be returned
  total.matching = empty.match.result(empty.rke())
  
  m = length(rke.list)  # no. of hospitals
  all.hids = rke.list.hospital.ids(rke.list)
  CHECK_EQ(all.hids, 1:m) # hospital ids should be 1, 2, 3...
  
  ## -1. New step: Match O-AB pairs first (i.e. ABO pairs)
  for(h in all.hids) {
    rke.h = rke.list[[h]]
    mOOU.h = OAB.matched.in.OUU(rke.h)
    # matchings are 3-way-only and of type OUU
    # make some sanity CHECKs here.
    CHECK_TRUE(mOOU.h$matching$utility %% 3 == 0)  # only OUU matches (3way)
    CHECK_MEMBER(mOOU.h$matching$match$pair.type, c("O", "U"))
    if(length(mOOU.h$OAB) > 0) {
      total.matching = add.matching(total.matching, mOOU.h$matching)
      rke.all <- rke.remove.pairs(rke.all, rm.pair.ids=c(mOOU.h$OAB, mOOU.h$U))
    }
  }
  
  ##  0. Compute IR constraints for S-subgraph.
  ir.constraints = compute.ir.constraints(pool, pair.types=c("S"), include.3way=T)
  
  # 1.  Match S internally -- Respect IR.
  s.subrke = rke.subgraph(rke.all, pair.type="S")
  s.constraints = subset(ir.constraints, pair.type=="S")
  s.matching = max.matching(s.subrke, ir.constraints=s.constraints,
                            include.3way=T)
  
  # 1a. Add matching information from S
  total.matching <- add.matching(total.matching, s.matching)
  ## remove some stuff that are not needed anymore
  rm(ir.constraints)
  
  ## 2. Matching R subgraph.
  ## 2a. Compute the overall and individual extended R-subgraphs
  #      and find 3-way max-matching internally.
  
  # AB, BA pairs (real and virtual)
  nAB = sapply(all.hids, function(h) nrow(subset(rke.all$pairs, desc=="A-B" & hospital==h)))
  nBA =  sapply(all.hids, function(h) nrow(subset(rke.all$pairs, desc=="B-A" & hospital==h)))
  virtual.nAB = rep(0, m)
  virtual.nBA = rep(0, m)
  # Compute the virtual AB/BA pairs.
  for(h in all.hids) {
    xRKE.h = rke.extended.Rsubgraph(rke.list[[h]])
    virtual.h = rke.count.virtual.pairs(xRKE.h)
    virtual.nAB[h] = virtual.h$AB
    virtual.nBA[h] = virtual.h$BA
  }
  # matches (ideal and real)
  mh.AB.ideal = sapply(all.hids, function(h) min(nAB[h], nBA[h] + virtual.nBA[h]))
  mh.BA.ideal = sapply(all.hids, function(h) min(nBA[h], nAB[h] + virtual.nAB[h]))

  # Unmatched R pairs (demand)
  zAB = nAB - mh.AB.ideal
  zBA = nBA - mh.BA.ideal
  # z_u (unmatched virtual pairs, kind of)
  uAB = sapply(all.hids, function(h) virtual.nAB[h] - max(0, mh.BA.ideal[h] - mh.AB.ideal[h]))
  uBA =sapply(all.hids, function(h) virtual.nBA[h] - max(0, mh.AB.ideal[h] - mh.BA.ideal[h]))
  
  # x_* = supply
  xAB = sum(zAB + uAB)
  xBA = sum(zBA + uBA)
  
  # y = allocation through the uniform share mechanism.
  yAB = rep(0, m)
  yBA = rep(0, m)
  if(sum(nAB) >= sum(nBA)) {
    yAB = g.share(zAB, xBA)
  } else {
    yBA = g.share(zBA, xAB)
  }
  # IR constraints
  xR.constraints = data.frame(pc=c(), pair.type=c(), hospital=c(), internal.matches=c())
  pcAB = subset(kPairs, desc=="A-B")$pc
  pcBA = subset(kPairs, desc=="B-A")$pc
  xRKE.all = rke.extended.Rsubgraph(rke.all)  # rke.all of the extended R subgraph
  
  # Iterate over hospitals
  # Calculates the mh.real  i.e. 
  # the internal matches of AB/BA pairs for each hospital h
  mh.AB.real = rep(0, m)  # real internal matches, m_h in text.
  mh.BA.real = rep(0, m)
  for(h in all.hids) {
    # get the extended R subgraph of hospital h
    xRKE.h = rke.extended.Rsubgraph(rke.list[[h]])
    # AB and BA pairs of hospital h
    real.R.pair.ids = subset(rke.all$pairs, hospital==h & pair.type=="R")$pair.id
    # compute max matching on the extended R subgraph of h
    m = max.matching(xRKE.h, promote.pair.ids=real.R.pair.ids, include.3way=T)
    
    mh.AB.real[h] = nrow(subset(m$match, desc=="A-B"))
    mh.BA.real[h] = nrow(subset(m$match, desc=="B-A"))
    xR.constraints = rbind(xR.constraints, 
                           list(pc=pcAB, pair.type="R", hospital=h, 
                                internal.matches=mh.AB.real[h]))
    xR.constraints = rbind(xR.constraints, 
                           list(pc=pcBA, pair.type="R", hospital=h, 
                                internal.matches=mh.BA.real[h]))
  }
  
  # deltas (difference between ideal and real)
  delta.AB = mh.AB.ideal - mh.AB.real
  delta.BA = mh.BA.ideal - mh.BA.real
  CHECK_TRUE(all(delta.AB >=0))
  CHECK_TRUE(all(delta.BA >=0))
  
  # 2d. Run the loop
  q = 0    
  loop.ended = F
  r.matching = empty.match.result(empty.rke())
  all.Rpair.ids = subset(rke.all$pairs, pair.type=="R")$pair.id
  while(!loop.ended) {
    loop.r.constraints = xR.constraints
    # in each loop (indexed by q) we get softer constraints.
    loop.r.constraints$internal.matches = 
      sapply(1:nrow(xR.constraints), 
             function(i) {
               hid = with(xR.constraints, hospital[i])
               pc = with(xR.constraints, pc[i])
               internal =  with(xR.constraints, internal.matches[i])
               if(pc==pcAB) {
                 return(internal + max(0, yAB[hid] + delta.AB[hid] - q))
               } else {
                 return(internal + max(0, yBA[hid] + delta.BA[hid] - q))
               }
             })
    r.matching = max.matching(xRKE.all, include.3way=T,
                              ir.constraints=loop.r.constraints,
                              promote.pair.ids=all.Rpair.ids)
    loop.ended = get.matching.status(r.matching) == "OK"
    q = q + 1
  }
  
  # 2e. Add matching information from R and remove pairs.
  total.matching <- add.matching(total.matching, r.matching)
  # rke.all = rke.remove.pairs(rke.all, rm.pair.ids=get.matching.ids(r.matching))
  
  #  3. Almost regular matching to the remainder
  ## Match OD's individually.
  hospital.ids <- rke.list.hospital.ids(rke.list)
  for(hid in hospital.ids) {
    rke.h = rke.list[[hid]]
    matched.hospital.ids <- intersect(get.matching.ids(total.matching),
                                      rke.pair.ids(rke.h))
    rke.h <- rke.remove.pairs(rke=rke.h,
                              rm.pair.ids=matched.hospital.ids)
    # promote UD pairs in this run.
    internal.matching = max.matching(rke.h, include.3way=T, regular.matching=T)
    # Add matching information from internal matching.
    total.matching <- add.matching(total.matching, internal.matching)
  }
  
  # match remainder.
  # We have to take the intersection because some matches (OUU) have been removed
  # from the entire graph (rke.all) -- see the top of xCM3
  rke.remainder = rke.remove.pairs(rke.all, rm.pair.ids=intersect(get.matching.ids(total.matching),
                                                                  rke.pair.ids(rke.all)))
  remainder.m = max.matching(rke.remainder, include.3way=T, regular.matching=T)
  total.matching <- add.matching(total.matching, remainder.m)
  
  CHECK_UNIQUE(total.matching$match$pair.id)
  return(total.matching)
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
  logdebug(sprintf("Running LOTTERY for UD=%d, theta=%d", ud.pc, theta))
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
  logdebug(sprintf("Before J-ball process. Sum=%d", get.sum()))
  logdebug("J bin")
  logdebug(set.J)
  while(length(set.J) > 0 && get.sum() < theta && plateau.count < 100) {
    h.sample = set.J[1]
    if(length(set.J) > 1)
      h.sample = sample(set.J, size=1)
    logdebug(sprintf(">> J ball: Hospital=%d", h.sample))
    # Draw without replacement.
    set.J = set.J[-which(set.J==h.sample)[1]]
    CHECK_MEMBER(Sh.XY[[h.sample]], Bh.XY[[h.sample]],
                 msg="Matched ids subset of all")
    
    avail = setdiff(Bh.XY[[h.sample]], Sh.XY[[h.sample]])
    if(length(avail) > 0) {
      xy.pair = uniform.sample(avail)
      logdebug(sprintf(">> J ball: Hospital=%d Adding pair %d (%s)", 
                       h.sample, xy.pair, subset(rke$pairs, pair.id==xy.pair)$desc))
      Sh.XY[[h.sample]] = c(Sh.XY[[h.sample]], xy.pair)
      plateau.count=0
    } else {
      plateau.count = plateau.count +1
      logdebug(sprintf(">> J ball: Hospital=%d -- NO PAIR AVAIL.", h.sample))
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
  logdebug("Computing MAX internal matchings..")
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
    rke.h = rke.list[[hid]]  # internal RKE
    # 2 . Compute a regular matching on the specific hospital (internal)
    internal.matching = max.matching(rke.h, include.3way=include.3way,
                                     regular.matching=F)
    ## Iterate over all UD codes
    for(upc in ud.pcs) {
      Qh[[hid]][upc] = nrow(subset(rke.h$pairs, pc==upc))
      CHECK_EQ(Qh[[hid]][[upc]], nrow(subset(rke.all$pairs, pc==upc & hospital==hid)))
      Sh[[hid]][[upc]] = subset(internal.matching$match, pc==upc)$pair.id
      logdebug(sprintf("Hospital %d is matching %d internally of %s", hid, length(Sh[[hid]][[upc]]),
                       subset(kPairs, pc==upc)$desc))
    }
  }
  ##  TO-DO(ptoulis): How do you know if you computed Qh, Sh correctly???
  return(list(Q=Qh, S=Sh))  
}

## Bonus mechanism. Ashlagi & Roth (2013)
Bonus = function(pool, include.3way=F) {
  CHECK_rke.pool(pool)
  # unload reports.
  rke.list = pool$rke.list
  rke.all = pool$rke.all
  
  # Total matching computed by the mechanism.
  total.matching = empty.match.result(empty.rke())
    
  ## 0. Initialize mechanism
  m = length(rke.list)  # no. of hospitals.
  CHECK_EQ(m %% 2, 0, msg="Run Bonus with even # of hospitals.")
  logdebug("######")
  logdebug("New BONUS")
  logdebug(sprintf("Total %d hospitals", m))
  
  ## 1. Match S pairs internally
  s.subrke = rke.subgraph(rke.all, pair.type="S")
  logdebug(sprintf("Computing max-matching in S-subgraph. Size %d pairs", rke.size(s.subrke)))
  s.matching = max.matching(s.subrke, include.3way=include.3way)
  
  logdebug("Performed S-matching..Utilities:")
  logdebug(get.matching.hospital.utilities(s.matching, m))
  
  ## 1b. Add S-matching
  total.matching <- add.matching(total.matching, s.matching)
  
  ## 2. Match R pairs. First match internally, then maximize global,
  #     promoting the internally-matched.
  promoted.Rpairs <- c()  # will contain the R pairs that were matched in internal M_R matchings.
  for(h in rke.list.hospital.ids(rke.list)) {
    rke.h = rke.list[[h]]
    # Maximize the matches on the R subgraph.
    h.Rpairs = subset(rke.h$pairs, pair.type=="R")$pair.id
    r.match.h <- max.matching(rke.h, include.3way=include.3way,
                              promote.pair.ids=h.Rpairs)
    matched.h.Rpairs = intersect(get.matching.ids(r.match.h), h.Rpairs)
    promoted.Rpairs <- c(promoted.Rpairs, matched.h.Rpairs)
  }
  r.subrke = rke.subgraph(rke.all, pair.type="R")
  r.matching = max.matching(r.subrke, include.3way=include.3way, 
                            promote.pair.ids=promoted.Rpairs)
  
  ## 2b. Add R-matching
  total.matching <- add.matching(total.matching, r.matching)
   
  #  3. Match OD/UD pairs. 
  #   3a. Split hospital in two sets
  k = as.integer(m/2)
  H.sets = list()
  H.sets[[1]] = sample(1:m, size=k, replace=F)
  H.sets[[2]] = setdiff(1:m, H.sets[[1]])

  logdebug("Split into two sets")
  logdebug(H.sets[[1]])
  logdebug(H.sets[[2]])
  
  # All Underdemanded PC codes
  ud.pcs = subset(kPairs, pair.type=="U")$pc
  logdebug("UD pair codes")
  logdebug(ud.pcs)
  
  # Compute Q-S:
  #   Q[hid][X-Y] = how many #X-Y in hospital hid
  #   S[hid][X-Y] = { ids of matched X-Y for hid }
  QS.obj = Bonus.QS(rke.pool=pool, include.3way=include.3way)
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
      logdebug("------------------")
      logdebug("NEW XY pair")
      logdebug("XY pairs from set=")
      logdebug(H.sets[[j]])
      tau.BHj.XY = subset(rke.all$pairs,
                          hospital %in% H.sets[[j]] & pair.id %in% all.XY)$pair.id
      # tau.BHother.XY = { ids that belong to H2 and are of type YX }
      tau.BHother.YX = subset(rke.all$pairs,
                              hospital %in% H.sets[[3-j]] & pair.id %in% all.YX)$pair.id
      
       # no. of such pairs (supply of YX pairs to be used in lottery)
      theta.j.YX = length(tau.BHother.YX)
      logdebug("X-Y pair=")
      logdebug(as.character(subset(kPairs, pc==i)$desc))
      logdebug("Y-X pair=")
      logdebug(as.character(subset(kPairs, pc==pc.reciprocal(i))$desc))
      logdebug("Before lottery: Pairs in J ball: ")
      x1 = sapply(H.sets[[j]], function(h) QS.obj$Q[[h]][i])
      x2 =sapply(H.sets[[j]], function(h) length(QS.obj$S[[h]][[i]]))
      x1.other = sapply(H.sets[[j]], function(hid) nrow(subset(rke.all$pairs, pc==i & hospital==hid)))
     
      logdebug("Printing (Q,S). Q=total #XY pairs,  S=total #XY that can be matched.")
      logdebug("Q=")
      logdebug(x1)
      logdebug("Some other way to get these numbers")
      logdebug(x1.other)
      logdebug("S=")
      logdebug(x2)
    
      logdebug("BEFORE lottery = ")
      logdebug("X-Y pairs set = {")
      logdebug(tau.BHj.XY)
      logdebug("Breakdown by hospital")
      logdebug(table(rke.pairs.hospitals(rke.all, pair.ids=tau.BHj.XY)))
      
      logdebug("Y-X pairs set = {")
      logdebug(tau.BHother.YX)
      logdebug("Breakdown by hospital")
      logdebug(table(rke.pairs.hospitals(rke.all, pair.ids=tau.BHother.YX)))
      
      
      ## Run the under-demanded lottery  (contains pair ids)
      S.XY = ud.lottery(rke=rke.all,
                        ud.pc=i,
                        Hn=H.sets[[j]],
                        theta=theta.j.YX,
                        QS=QS.obj)
    
      
      logdebug("AFTER lottery = ")
      logdebug("X-Y pairs matched in lottery = {")
      logdebug(S.XY)
      logdebug("Breakdown by hospital")
      logdebug(table(rke.pairs.hospitals(rke.all, pair.ids=S.XY)))
      
      CHECK_DISJOINT(S.XY, tau.BHother.YX, msg="Not even the same type")
      CHECK_MEMBER(S.XY, tau.BHj.XY, msg="Matched subset of all")
      ## the only ids to be considered in the matching
      XYYX.ids = union(S.XY, tau.BHother.YX)
      # Here you match   X-Y pairs from Hj   with   Y-X pairs from Hother
      subrke = rke.keep.pairs(rke.all, pair.ids=XYYX.ids)
      ## This is an interesting point.
      Mj = max.matching(subrke, include.3way=include.3way, regular.matching=F)
      
      logdebug("Performed OU-matching..Utilities:")
      logdebug(get.matching.hospital.utilities(Mj, m))
      
      ## 0.3 Add new matching
      total.matching <- add.matching(total.matching, Mj)
    }
  }
  return(total.matching)
}