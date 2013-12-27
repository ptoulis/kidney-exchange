##  Tables of the paper. 
## Code to save into latex table.
library(xtable)

rm(list=ls())
# Add necessary libs.
source("../r-toolkit/checks.R")
source("../r-toolkit/logs.R")
source("terminology.R")
source("rke.R")
source("matching.R")
source("mechanisms.R")

## hand-made bootstrap (USDA organic)
bootstrap.mean = function(x) sd(replicate(1000, { y = sample(x, replace=T); mean(y)}))

## Combines the result matrix with the SE errors.
add.se = function(M, SE) {
    D = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    colnames(D) = colnames(M)
    for(i in 1:nrow(M)) {
        D[i,]= sapply(1:ncol(M), function(j) sprintf("%.2f (%.2f)", M[i,j], SE[i,j]))
    }
    return(D)
}

create.comparison <- function(mechanisms,
                              nHospitals, nSize,
                              uniform.pra, include.3way,
                              baseline.strategy,
                              deviation.strategy,
                              nsamples) {
  # Comparison is an object that defines what we are to compare.
  # Given the arguments the comparison defines:
  #  1) Which mechanisms to compare, i.e. mechanisms=c("rCM", "xCM", etc)
  #  2) Hospital setup, i.e. #hospitals, #pairs/hospitals
  #  3) Compatibility/matching setup, i.e. uniform PRA, and include or not 3way
  #  4) Strategies to compare, i.e. baseline="tttt...", deviation="ctttt..."
  #  5) #trials to take
  CHECK_MEMBER(mechanisms, c("rCM", "xCM", "Bonus"))
  CHECK_GE(nHospitals, 1)
  CHECK_GE(nSize, 1)
  CHECK_GE(nsamples, 1)
  if(nsamples < 10)
    warning("Maybe too few samples?")
  baseline.str = strsplit(baseline.strategy, split="")[[1]]
  deviation.str = strsplit(deviation.strategy, split="")[[1]]
  CHECK_MEMBER(baseline.str, c("t", "c", "r"))
  CHECK_MEMBER(deviation.str, c("t", "c", "r"))
  return(list(mechanisms=mechanisms, 
              m=nHospitals, n=nSize, 
              uniform.pra=uniform.pra, include.3way=include.3way,
              baseline.strategy=baseline.strategy,
              deviation.strategy=deviation.strategy,
              nsamples=nsamples))
}

CHECK_comparison <- function(comparison) {
  CHECK_MEMBER(c("baseline.strategy", "deviation.strategy", "m", "n", 
                 "uniform.pra", "include.3way", "mechanisms"),
               names(comparison))
  CHECK_MEMBER(comparison$include.3way, c(T, F))
  CHECK_MEMBER(comparison$uniform.pra, c(T, F))
  # Check whether tttt == tttc ... in length
  CHECK_EQ(nchar(comparison$baseline.strategy), nchar(comparison$deviation.strategy))
  # Check whether the strategy profile string has the correct length.
  CHECK_EQ(nchar(comparison$baseline.strategy), comparison$m)
}

# An example COMPARISON object.
example.comparison = create.comparison(mechanisms=c("rCM", "xCM", "Bonus"),
                                       nHospitals=4, nSize=35,
                                       uniform.pra=T, include.3way=F,
                                       baseline.strategy="tttt",
                                       deviation.strategy="cccc",
                                       nsamples=100)

compare.mechanisms <- function(comparison) {
  # Compares two mechanisms.
  # Result looks like  results[[strategy]][[mechanism]][[..attr..]]
  #
  # attr= { utility, hospital.pairsBreakdown, 
  #         total.matchingInfo, internal.matchingInfo, mech.matchingInfo}
  # For example,
  #   result$baseline$rCM$utility = (m x samples) utility matrix
  #   result$baseline$rCM$internal.matchingInfo[[3]] = MATCHING_INFO for hospital 3
  #   result$baseline$rCM$mech.matchingInfo = MATCHING INFO for only the mechanism
  CHECK_comparison(comparison)
  all.strategies <- c("baseline", "deviation")
  result = list(baseline=list(strategy=comparison$baseline.strategy),
                deviation=list(strategy=comparison$deviation.strategy))
  
  empty.matchingInfo <- empty.match.result(empty.rke())$information
  empty.pairsBreakdown <- table(subset(kPairs, pc < 0)$pair.type)
  # 0. Initialize the result object.
  for (mech in comparison$mechanisms) {
    for (s in all.strategies) {
      # Setting the result object for (strategy, mechanism)
      result[[s]][[mech]] <- list(utility=matrix(0, nrow=comparison$m, ncol=comparison$nsamples),
                                  mech.matchingInfo=empty.matchingInfo,
                                  internal.matchingInfo=list(),
                                  hospital.pairsBreakdown=list(),
                                  total.matchingInfo=empty.matchingInfo)
      
      for (hid in 1:comparison$m) {
        # Initialize the matching info
        result[[s]][[mech]]$internal.matchingInfo[[hid]] = empty.matchingInfo
        # set empty table of pair types (O, R, S, U)
        result[[s]][[mech]]$hospital.pairsBreakdown[[hid]] = empty.pairsBreakdown
      }
    }
  }  # done with initialization
  
  pb = txtProgressBar(style=3)
  for(i in 1:comparison$nsamples) {
    # 1. Sample the RKE pool
    rke.pool = rrke.pool(m=comparison$m, 
                         n=comparison$n,
                         uniform.pra=comparison$uniform.pra)
    # 2. Create the KPD markets
    kpds = list()
    kpds$baseline = kpd.create(rke.pool=rke.pool,
                               strategy.str=comparison$baseline.strategy)
    kpds$deviation = kpd.create(rke.pool=rke.pool,
                                strategy.str=comparison$deviation.strategy)
    
    infoMap = list(total.matchingInfo="total.matching",
                   mech.matchingInfo="mech.matching")
    # 3. Main loop starts here.
    for (mech in comparison$mechanisms) {
      for(s in all.strategies) {
        # 3.1 Run mechanism -> matching output for a specific strategy profile
        match =  Run.Mechanism(kpd=kpds[[s]], mech=mech, include.3way=comparison$include.3way) 
        # this object has (total.matching, mech.matching, internal.matchings)
        
        # 3.2 Update utility matrix.
        result[[s]][[mech]]$utility[, i] <- get.matching.hospital.utilities(match$total.matching,
                                                                            comparison$m)
        # 3.3 Update total + mech information
        for(info in names(infoMap)) {
          mappedTo = infoMap[[info]]
          result[[s]][[mech]][[info]] = result[[s]][[mech]][[info]] + match[[mappedTo]]$information
        }  # for info types
        # 3.4 Update hospital internal matching info
        for (hid in 1:comparison$m) {
          result[[s]][[mech]]$internal.matchingInfo[[hid]] =  
            result[[s]][[mech]]$internal.matchingInfo[[hid]] + match$internal.matchings[[hid]]$information
          ## pairs breakdown
          result[[s]][[mech]]$hospital.pairsBreakdown[[hid]] = 
            result[[s]][[mech]]$hospital.pairsBreakdown[[hid]] + 
            table(subset(match$total.matching$match, hospital==hid)$pair.type)
        }  # for all hospitals, matching info
      } # all strategies (baseline, deviation)
      
    } # for every mechanism
    setTxtProgressBar(pb, value=i / comparison$nsamples)
  }  # all trials
  
  ## Return the final result
  return(result)                       
}

summarize.compare.output <- function(compare.mech.out) {
  ## Gives a verbose summary of the outcome of the comparison
  #
  # Args:
  #   compare.mech.out = output of compare.mechanisms()
  all.mechanisms = c("rCM", "xCM", "Bonus")
  for(strategy in names(compare.mech.out)) {
    for(mech in intersect(names(compare.mech.out[[strategy]]), all.mechanisms)) {
      print(sprintf("Strategy type= %s, strategy profile= %s, mechanism= %s, ", 
                    strategy, compare.mech.out[[strategy]]$strategy, mech))
      print("Hospital utilities (avg.)")
      print(rowMeans(compare.mech.out[[strategy]][[mech]]$utility))
      print("Total welfare (sum)")
      print(sum(compare.mech.out[[strategy]][[mech]]$utility))
    }
  }
}

mech.weakness <- function(mechanism, nHospitals=6, nSize=25, ntrials=10) {
  CHECK_EQ(nHospitals %% 2, 0)
  kCurrentLogLevel <<- 5
  OU = which(kPairs$desc %in% c("O-A", "A-O"))
  kPairs$prob[-OU] <<- 0
  print(kPairs)
  ud.pairs = nrow(subset(rrke(1000)$pairs, pair.type=="U"))
  print(binom.test(x=ud.pairs, n=1000, p=5/6))
  ##
  mech=mechanism
  utils.t <- c()
  utils.c <- c()
  pb <- txtProgressBar(style=3)
  
  for(i in 1:ntrials) {
    pool = rrke.pool(m=nHospitals, n=nSize, uniform.pra=T)
    str.truth <- paste(rep("t", nHospitals), collapse="")
    str.dev <- paste(c("c", rep("t", nHospitals-1)), collapse="")
    kpd.t <- kpd.create(pool, strategy.str=str.truth)
    kpd.c <- kpd.create(pool, strategy.str=str.dev)
    
    kLogFile <<- "Bonus-H1-truthful.log"
    xt = Run.Mechanism(kpd.t, mech=mech,  include.3way=F)
    utils.t[i] <- get.matching.hospital.utilities(xt$total.matching, nHospitals)[1]
    kLogFile <<- "Bonus-H1-deviates.log"
    xc = Run.Mechanism(kpd.c, mech=mech,  include.3way=F)
    utils.c[i] <- get.matching.hospital.utilities(xc$total.matching, nHospitals)[1]
    
    setTxtProgressBar(pb, value=i/ntrials)
    print(sprintf("Truthful H1 (%s)", mech))
    print(summary(utils.t))
    print(sprintf("Deviating H1 (%s)", mech))
    print(summary(utils.c))
  }
  
  print(summary(utils.t))
  print("Deviation")
  print(summary(utils.c))
}

mech.weakness.theoretical <- function(nHospitals=6, nSize=25, ntrials=10) {
  CHECK_EQ(nHospitals %% 2, 0)
  util.t <- c()
  util.c <- c()
  bad.cases <- 0
  for(i in 1:ntrials) {
    U = rbinom(nHospitals, size=nSize, prob=5/6)
    O = nSize-U
    theta = sum(tail(O, nHospitals/2))  # supply of over-demanded.
    if(any(U < O)) {
      bad.cases <- bad.cases + 1
      next
    }
    # H1=truthful case.
    m.XY.t <- head(O, nHospitals/2)  # matched XY from side 1 = #YX pairs (PM assumption)
    S1.t = sum(m.XY.t)  # XY pairs of one side that are matched
    reports.XY.t <- head(U, nHospitals/2)
    excess.t = theta - S1.t  # excess of supply to be allocated
    excess.H1.t =  excess.t * reports.XY.t[1] / sum(reports.XY.t)  # excess that could go to H1
    
    # deviation case.
    m.XY.c <- head(O, nHospitals/2) # matched XY from side 1
    m.XY.c[1] <- 0
    S1.c = sum(m.XY.c)  # XY pairs of one side that are matched
    reports.XY.c <- head(U, nHospitals/2)
    reports.XY.c[1] <- U[1] - O[1]
    excess.c = theta - S1.c  # excess of supply to be allocated
    excess.H1.c =  excess.c * reports.XY.c[1] / sum(reports.XY.c)  # excess that could go to H1
    
    # print(excess.H1.t)
    # print(excess.H1.c)
    util.t <- c(util.t, O[1] + max(excess.H1.t, 0) + O[1])
    util.c = c(util.c, max(excess.H1.c, 0) + 2 * O[1])    
  }
  print(sprintf("Bad cases = %d / %d", bad.cases, ntrials))
  print("Truthful case...")
  print(summary(util.t))
  print("Deviation case..")
  print(summary(util.c))
}

table.regularity <- function(nsamples, include.3way, verbose=F) {
  # Explore regularity assumptions for 2-way or 3-way exchanges
  #
  # Regularity has "aspects" which are basically things we check to see
  # whether regularity holds.
  # The function count.aspect(rke, aspect) is the basic function here
  # and gives the value for this "aspect" for a particular RKE object.
  # This experiment computes and stores aspect values for the original RKE
  # and the remainder RKE after a structured matching (e.g. regular has been imposed)
  #
  # For 2-way exchanges:
  #   aspect R = count #short side of R (bc in theory short side=0 after regular matching)
  #   aspect S = count #S pairs (bc. S pairs=0,1 after regular matching)
  #   aspect O = count # pairs  (bc. O pairs=0 after regular matching)
  #
  # For 3-way exchanges
  #   aspect R = sum of R pairs  (bc. R pairs=0 after extended R max matching)
  #   aspect S= sum of S pairs (b.c. S pairs =0 after 3way matching)
  #   aspect O = 0
  #
  all.sizes <- round(seq(5, 200, by=25))
  if(include.3way)
    all.sizes = round(seq(10, 140, by=20))  # smaller sizes for 3-way
  
  sampled.sizes <- sample(all.sizes, size=nsamples, replace=T)
  print("Sampled sizes breakdown")
  print(table(sampled.sizes))
  
  # stores information about the matches e.g. OR, OO, OS, ...
  regularity.aspects = c("nsize", "O", "R", "S", "O.unmatched", "R.unmatched", "S.unmatched")
  regularity = list(UPRA=matrix(NA, nrow=0, ncol=length(regularity.aspects)),
                    NonUPRA=matrix(NA, nrow=0, ncol=length(regularity.aspects)))
  colnames(regularity$UPRA) = regularity.aspects
  colnames(regularity$NonUPRA) = regularity.aspects
  
  count.aspect <- function(rke, aspect) {
    # Counts the aspect of the regularity assumption 
    # we are interested in.
    # If aspect=R then count the size of the short R side
    # If aspect=O count #O
    # If aspect=S count S's (actually max(0, n-1)
    CHECK_MEMBER(aspect, c("O", "R", "S"))
    if(aspect=="R") {
      nAB = nrow(subset(rke$pairs, desc=="A-B"))
      nBA = nrow(subset(rke$pairs, desc=="B-A"))
      if(include.3way) {
        # for 3-way we need to check whether all R pairs were matched
        return(nAB+nBA)
      } else {
        # for 2-way we need to check whether the short side is matched.
        return(min(nAB, nBA))
      }
    } else if(aspect=="O") {
      # for 3-way we don't care about O pairs.
      if(include.3way) return(0)
      # for 2-way we check the #O pairs unmatched.
      # same story for S pairs.
      return(nrow(subset(rke$pairs, pair.type=="O")))
    } else {
      if(include.3way) {
        return(nrow(subset(rke$pairs, pair.type=="S")))
      }
      stypes = c("O-O", "A-A", "B-B", "AB-AB")
      return(sum(sapply(stypes, function(t) {
        ns = nrow(subset(rke$pairs, desc==t))
        return(max(0, ns-1))
      })))
    }
  }
  
  pb = txtProgressBar(style=3)
  CHECK_EQ(length(sampled.sizes), nsamples)
  
  # save at the last iteration too.
  save.checkpoints <- c(seq(1, nsamples-1, length.out=20), nsamples)
  for(uniform.pra in c(T, F)) {
    print("")
    print(sprintf("PRA=%s, nsamples=%d, 3way=%s", uniform.pra, nsamples, include.3way))
    for(i in 1:nsamples) {
      n = sampled.sizes[i]
      rke = rrke(n, uniform.pra=uniform.pra)
      ## compute the regular 
      Rpair.ids = rke.filter.pairs(rke, attr="pair.type", value="R")
      Spair.ids = rke.filter.pairs(rke, attr="pair.type", value="S")
      CHECK_DISJOINT(Rpair.ids, Spair.ids)
      OU.pair.ids = setdiff(rke.pair.ids(rke), union(Rpair.ids, Spair.ids))
      
      update <- NA  # this will be the row to be added to results
      if(include.3way) {
        xRKE = rke.extended.Rsubgraph(rke)      
        m3way = max.matching(xRKE, include.3way=T, promote.pair.ids=Rpair.ids)
        # computes unmatched AB + BA
        unmatched.R = count.aspect(rke.remove.pairs(xRKE, get.matching.ids(m3way)), "R")
        aspect.R = 0
        nAB = nrow(subset(xRKE$pairs, desc=="A-B"))
        nBA = nrow(subset(xRKE$pairs, desc=="B-A"))
        CHECK_EQ(nAB + nBA, length(Rpair.ids))
        nmatches.real = nAB + nBA - unmatched.R  # matches realized
        excess = abs(nBA-nAB)
        short.side = min(nAB, nBA)
        virtualR = rke.count.virtual.pairs(xRKE)
        
        nmatches.under.regularity = 0
        
        if(nAB > nBA) {
          nmatches.under.regularity = 2 * short.side + min(excess, virtualR$BA)
        } else {
          nmatches.under.regularity = 2 * short.side + min(excess, virtualR$AB)
        }
        # This will ne included in the violations row update.
        aspect.R = nmatches.under.regularity - nmatches.real
        if(verbose) {
          print("")
          plot(xRKE)
          print("Max matching on extended R subgraph")
          print(m3way$matched.cycles)

          print(sprintf("Total AB=%d BA=%d Virtual AB=%d, BA=%d, unmatched R pairs=%d  Violations=%d",
                        nAB, nBA, virtualR$AB, virtualR$BA, unmatched.R, aspect.R))
          readline("Press ENTER")
        }
        
        # S-subgraph violations.
        sRKE = rke.keep.pairs(rke, pair.ids=Spair.ids)
        mS = max.matching(sRKE, include.3way=T)
        update = c(n, 
                   0,
                   count.aspect(rke, "R"),
                   count.aspect(rke, "S"),
                   0,
                   aspect.R,
                   count.aspect(rke.remove.pairs(rke, get.matching.ids(mS)), "S"))
      } else {
        OU.rke = rke.keep.pairs(rke, pair.ids=OU.pair.ids)
        S.rke = rke.keep.pairs(rke, pair.ids=Spair.ids)
        R.rke = rke.keep.pairs(rke, pair.ids=Rpair.ids)
        mOU = max.matching(OU.rke)
        mR = max.matching(R.rke)
        mS = max.matching(S.rke)
        if(verbose) {
          print("")
          print(sprintf("There are %d AB %d BA, Matched %d AB and %d BA", 
                        nrow(subset(rke$pairs, desc=="A-B")),
                        nrow(subset(rke$pairs, desc=="B-A")),
                        nrow(subset(mR$match, desc=="A-B")),
                        nrow(subset(mR$match, desc=="B-A"))))
        }
        # Update regularity matrix
        update = c(n, 
                   count.aspect(rke, "O"),
                   count.aspect(rke, "R"),
                   count.aspect(rke, "S"),
                   count.aspect(rke.remove.pairs(OU.rke, get.matching.ids(mOU)), "O"),
                   count.aspect(rke.remove.pairs(R.rke, get.matching.ids(mR)), "R"),
                   count.aspect(rke.remove.pairs(S.rke, get.matching.ids(mS)), "S"))
      }
      
      if(uniform.pra) {
        regularity$UPRA = rbind(regularity$UPRA, update)
        rownames(regularity$UPRA) <- NULL
      } else {
        regularity$NonUPRA = rbind(regularity$NonUPRA, update)
        rownames(regularity$NonUPRA) <- NULL
      }
      
      setTxtProgressBar(pb, value=i/nsamples)
      filename = ""
      if(i %in% save.checkpoints) {
        if(include.3way) {
          filename = "out/table1-3way.Rdata"
          table1.3way = regularity
          save(table1.3way, file=filename)
        } else {
          filename = "out/table1-2way.Rdata"
          table1.2way = regularity
          save(table1.2way, file=filename)
        }
      }
    }
  }
  print("")
  print(sprintf("Simulation complete. File saved in %s", filename))
}

table.welfare.incentives <- function(nhospitals=6, nsize=15, 
                                     include.3way=F, nsamples=100) {
  # Table of 2way exchanges to compare Welfare and Incentives.
  # or Table of 3way exchanges to compare welfare + incentives
  results <<- list()
  get.strategy.profile <- function(no.truthful) {
    CHECK_TRUE(no.truthful >= 0 & no.truthful <= nhospitals, msg="#truthful should be correct")
    no.deviating = nhospitals - no.truthful
    return(paste(c(rep("t", no.truthful), rep("c", no.deviating)), collapse=""))
  }
  
  run.comparison <- function(base.Nt, dev.Nt, pra) {
    # Runs a single experiment.
    # 
    # Args:
    #  base.Nt = no. of truthful agents in baseline strategy profile
    #  dev.Nt = no. of truthful strategies in dev strategy profile
    #  pra = T or F, whether we want uniform PRA or non-uniform PRA.
    baseline.strategy = get.strategy.profile(base.Nt)
    deviation.strategy = get.strategy.profile(dev.Nt)
    print("")
    print(sprintf("Comparing profiles  %s vs. %s, PRA=%s, 3way=%d", 
                  baseline.strategy, deviation.strategy, pra,
                  include.3way))
    comparison = create.comparison(mechanisms=c("rCM", "xCM", "Bonus"),
                                   nHospitals=nhospitals, nSize=nsize,
                                   uniform.pra=pra, 
                                   include.3way=include.3way,
                                   baseline.strategy=baseline.strategy,
                                   deviation.strategy=deviation.strategy,
                                   nsamples=nsamples)
    profile.name = sprintf("prof%d%d", base.Nt, dev.Nt)
    result.name = sprintf("%s.%s.%s", 
                          profile.name,
                          ifelse(pra, "UPRA", "NonUPRA"),
                          ifelse(include.3way, "3way", "2way"))
    results[[result.name]] <<- compare.mechanisms(comparison)
    save(results, file=sprintf("out/table%d-m%dn%d-results.Rdata", 
                               ifelse(include.3way, 3, 2),
                               nhospitals, nsize))
  }
  
  run.comparison(nhospitals, nhospitals-1, T)
  run.comparison(nhospitals, nhospitals-1, F)
  run.comparison(1, 0, T)
  run.comparison(1, 0, F)  
}