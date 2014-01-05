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
  CHECK_MEMBER(mechanisms, kImplementedKPDMechanisms)
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
                                       nHospitals=4, nSize=18,
                                       uniform.pra=T, include.3way=F,
                                       baseline.strategy="tttt",
                                       deviation.strategy="cccc",
                                       nsamples=10)

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
  # main loop
  for(i in 1:comparison$nsamples) {
    # 1. Sample the RKE pool
    rke.pool = rrke.pool(m=comparison$m, 
                         n=comparison$n,
                         uniform.pra=comparison$uniform.pra)
    # 2. Create the KPD markets
    kpds = list()
    kpds$baseline = kpd.create(rke.pool=rke.pool,
                               strategy.str=comparison$baseline.strategy,
                               include.3way=comparison$include.3way)
    kpds$deviation = kpd.create(rke.pool=rke.pool,
                                strategy.str=comparison$deviation.strategy,
                                include.3way=comparison$include.3way)
    
    # 3. Main loop starts here.
    for (mech in comparison$mechanisms) {
      for(s in all.strategies) {
        # 3.1 Run mechanism -> matching output for a specific strategy profile
        run.output =  Run.Mechanism(kpd=kpds[[s]], mech=mech, include.3way=comparison$include.3way) 
        # this object has (total.matching, mech.matching, internal.matchings)
        
        # 3.2 Update utility matrix.
        result[[s]][[mech]]$utility[, i] <- get.matching.hospital.utilities(run.output$total.matching,
                                                                            comparison$m)
        # 3.3 Update total + mech information
        for(m in c("total.matching", "mech.matching")) {
          key = sprintf("%sInfo", m)
          result[[s]][[mech]][[key]] = result[[s]][[mech]][[key]] + run.output[[m]]$information
        }
        # 3.4 Update hospital internal matching info
        for (hid in 1:comparison$m) {
          result[[s]][[mech]]$internal.matchingInfo[[hid]] =  
            result[[s]][[mech]]$internal.matchingInfo[[hid]] + run.output$internal.matchings[[hid]]$information
          ## pairs breakdown
          result[[s]][[mech]]$hospital.pairsBreakdown[[hid]] = 
            result[[s]][[mech]]$hospital.pairsBreakdown[[hid]] + 
            table(subset(run.output$total.matching$match, hospital==hid)$pair.type)
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

safe.sample <- function(from, size) {
  CHECK_TRUE(length(from) >= size)
  ind = sample(c(rep(1, size), rep(0, length(from)-size)))
  return(from[which(ind==1)])
}

table.regularity <- function(nsamples, max.hospitalSize=500,
                             include.3way, verbose=F) {
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
  all.sizes <- round(seq(10, max.hospitalSize, by=50))
  if(include.3way)
    all.sizes = round(seq(10, max.hospitalSize, by=50))  # smaller sizes for 3-way
  
  sampled.sizes <- sample(all.sizes, size=nsamples, replace=T)
  print("Sampled sizes breakdown")
  print(table(sampled.sizes))
  
  # stores information about the matches e.g. OR, OO, OS, ...
  # regularity.aspects = c("nsize", "O", "R", "S", "O.unmatched", "R.unmatched", "S.unmatched")
  S.aspects = c("A-A", "B-B", "O-O", "AB-AB")
  R.aspects = c("R")
  OU.aspects = c("O")
  if(include.3way) OU.aspects = c("O", "O-AB")
  
  unmatched.aspects <- function(x) {
    as.vector(sapply(x, function(s) sprintf("%s.unmatched", s)))  
  }

  regularity.aspects = c("nsize",
                         S.aspects, unmatched.aspects(S.aspects),
                         R.aspects, unmatched.aspects(R.aspects),
                         OU.aspects, unmatched.aspects(OU.aspects))
  regularity = list(UPRA=matrix(NA, nrow=0, ncol=length(regularity.aspects)),
                    NonUPRA=matrix(NA, nrow=0, ncol=length(regularity.aspects)))
  colnames(regularity$UPRA) = regularity.aspects
  colnames(regularity$NonUPRA) = regularity.aspects

  print(regularity.aspects)
  
  count.aspects <- function(rke, aspects) {
    sapply(aspects, function(what) {
      if(what %in% kPairs$pair.type) {
        # if count a specific pair type, e.g. S or R
        return(nrow(subset(rke$pairs, pair.type==what)))
      } else if(what %in% kPairs$desc) {
        # if count a specific donor-patient type (e.g. A-B)
        return(nrow(subset(rke$pairs, desc==what)))
      } else {
        stop(sprintf("Could not find aspect %s", what))
      }
    })
  }
  
  pb = txtProgressBar(style=3)
  CHECK_EQ(length(sampled.sizes), nsamples)
  
  # save at the last iteration too.
  save.checkpoints <- c(seq(1, nsamples-1, length.out=20), nsamples)
  for(uniform.pra in c(T, F)) {
    print("")
    print(sprintf("PRA=%s, nsamples=%d, 3way=%s", uniform.pra, nsamples, include.3way))
    for(i in 1:nsamples) {
      # size of the RKE object
      nsize = sampled.sizes[i]
      rke = rrke(nsize, uniform.pra=uniform.pra)
      ## Pair ids (common to both 2-way and 3-way loops)
      Rpair.ids = rke.filter.pairs(rke, attr="pair.type", value="R")
      ABpair.ids = subset(rke$pairs, desc=="A-B")$pair.id
      BApair.ids = subset(rke$pairs, desc=="B-A")$pair.id
      Spair.ids = rke.filter.pairs(rke, attr="pair.type", value="S")
      OUpair.ids = setdiff(rke.pair.ids(rke), union(Rpair.ids, Spair.ids))
      CHECK_DISJOINT(Rpair.ids, Spair.ids)
      
      update <- NA  # this will be the row to be added to results
      only.2way = 1-include.3way
      
      if(only.2way) {
        # 2-way exchanges
        # Preprocess S-rke (updates Spair.ids)
        # For all T-T remove one pair at random if it is odd.
        for(aspect in S.aspects) {
          TTpair.ids = subset(rke$pairs, desc==aspect)$pair.id
          if(length(TTpair.ids) > 0 && length(TTpair.ids) %% 2 == 1) {
            rm.id = safe.sample(TTpair.ids, size=1)
            if(verbose) {
              print(sprintf("Graph %s had %d pairs. Removing one.",
                            aspect, length(TTpair.ids)))
            }
            Spair.ids = setdiff(Spair.ids, rm.id)
          }
        }
        #
        # Pre-process R subgraph (updates Rpair.ids)
        # First balance-then match.
        excess = abs(length(ABpair.ids) - length(BApair.ids))
        # remove excess.
        if(length(ABpair.ids)  > length(BApair.ids)) {
          Rpair.ids = setdiff(Rpair.ids, safe.sample(ABpair.ids, size=excess))
        } else {
          Rpair.ids = setdiff(Rpair.ids, safe.sample(BApair.ids, size=excess))
        }
        # No pre-process for OU matches.
        OU.rke = rke.keep.pairs(rke, pair.ids=OUpair.ids)
        S.rke = rke.keep.pairs(rke, pair.ids=Spair.ids)
        R.rke = rke.keep.pairs(rke, pair.ids=Rpair.ids)
        
        mS = max.matching(S.rke)
        S.rke.remainder = rke.remove.pairs(S.rke, get.matching.ids(mS))
        mR = max.matching(R.rke)
        R.rke.remainder = rke.remove.pairs(R.rke, get.matching.ids(mR))
        mOU = max.matching(OU.rke)
        OU.rke.remainder = rke.remove.pairs(OU.rke, get.matching.ids(mOU))
              
        # Update regularity matrix
        update = as.vector(c(nsize, 
                   count.aspects(S.rke, S.aspects),
                   count.aspects(S.rke.remainder, S.aspects), 
                   count.aspects(R.rke, R.aspects),
                   count.aspects(R.rke.remainder, R.aspects),
                   count.aspects(OU.rke, OU.aspects),
                   count.aspects(OU.rke.remainder, OU.aspects)))
      }
      else if(include.3way) {
        ## No pre-processing for S-subgraph
        
        # Preprocess for R-subgraph
        xRKE = rke.extended.Rsubgraph(rke)   
        virtualR = rke.count.virtual.pairs(xRKE)
        if(length(ABpair.ids) < length(BApair.ids)) {
          if(length(ABpair.ids) + virtualR$AB < length(BApair.ids)) {
            rm.ids = safe.sample(from=BApair.ids, size=length(BApair.ids) - length(ABpair.ids) - virtualR$AB)
            Rpair.ids = setdiff(Rpair.ids, rm.ids)
          }
        } else {
          if(length(BApair.ids) + virtualR$BA < length(ABpair.ids)) {
            rm.ids = safe.sample(from=ABpair.ids,
                                 size=length(ABpair.ids) - length(BApair.ids) - virtualR$BA)
            Rpair.ids = setdiff(Rpair.ids, rm.ids)
          }
        }
        # Pre-process for OU graph
        # See if OUU match all O-AB (i.e. AB-O with new notation)
        mOUU = OAB.matched.in.OUU(rke)
        virtualR = rke.count.virtual.pairs(xRKE, return.bipartite=T)
        virtual.pair.ids = c(virtualR$AB, virtualR$BA)
        select.rm.pair.ids = rbinom(length(virtual.pair.ids), size=1, prob=0.5)
        rm.pair.ids = virtual.pair.ids[which(select.rm.pair.ids==1)]
        # Remove the OUU matches and the random virtual pairs.
        OUpair.ids = setdiff(OUpair.ids, c(rm.pair.ids, mOUU$U, mOUU$OAB))
        
        # Compute the RKE after the transformations.
        S.rke = rke.keep.pairs(rke, Spair.ids)
        R.rke = xRKE
        OU.rke = rke.keep.pairs(rke, OUpair.ids)
        
        mS = max.matching(S.rke, include.3way=T)
        S.rke.remainder = rke.remove.pairs(S.rke, get.matching.ids(mS))
        mR = max.matching(R.rke, include.3way=T, promote.pair.ids=Rpair.ids)
        R.rke.remainder = rke.remove.pairs(R.rke, get.matching.ids(mR))
        mOU = max.matching(OU.rke, regular.matching=T)  # should be 2-way
        OU.rke.remainder = rke.remove.pairs(OU.rke, get.matching.ids(mOU))
        
        # Update regularity matrix
        update = as.vector(c(nsize, 
                             count.aspects(S.rke, S.aspects),
                             count.aspects(S.rke.remainder, S.aspects), 
                             count.aspects(R.rke, R.aspects),
                             count.aspects(R.rke.remainder, R.aspects),
                             count.aspects(OU.rke, OU.aspects),
                             count.aspects(OU.rke.remainder, OU.aspects)))
      } else {
        # this is a 2-way exchange
        stop("No include.3way?")
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


table.welfare.incentives <- function(mechanisms=kImplementedKPDMechanisms,
                                     nhospitals=6, nsize=15, 
                                     include.3way=F, nsamples=100,
                                     filename.prefix="table",
                                     run.profiles=c(1,2,3,4)) {
  # Table of 2way exchanges to compare Welfare and Incentives.
  # or Table of 3way exchanges to compare welfare + incentives

  get.strategy.profile <- function(no.truthful) {
    CHECK_TRUE(no.truthful >= 0 & no.truthful <= nhospitals, msg="#truthful should be correct")
    no.deviating = nhospitals - no.truthful
    return(paste(c(rep("t", no.truthful), rep("c", no.deviating)), collapse=""))
  }
  
  run.comparison <- function(results, base.Nt, dev.Nt, pra) {
    # Runs a single experiment.
    # 
    # Args:
    #  base.Nt = no. of truthful agents in baseline strategy profile
    #  dev.Nt = no. of truthful strategies in dev strategy profile
    #  pra = T or F, whether we want uniform PRA or non-uniform PRA.
    baseline.strategy = get.strategy.profile(base.Nt)
    deviation.strategy = get.strategy.profile(dev.Nt)
    print("")

    comparison = create.comparison(mechanisms=mechanisms,
                                   nHospitals=nhospitals, nSize=nsize,
                                   uniform.pra=pra, 
                                   include.3way=include.3way,
                                   baseline.strategy=baseline.strategy,
                                   deviation.strategy=deviation.strategy,
                                   nsamples=nsamples)
    
    print(sprintf("Comparing profiles  %s vs. %s, PRA=%s, 3way=%d, m=%d, n=%d", 
                  baseline.strategy, deviation.strategy, pra,
                  include.3way,
                  comparison$m,
                  comparison$n))
    
    profile.name = sprintf("prof%d%d", base.Nt, dev.Nt)
    result.name = sprintf("%s.%s.%s", 
                          profile.name,
                          ifelse(pra, "UPRA", "NonUPRA"),
                          ifelse(include.3way, "3way", "2way"))
    results[[result.name]] <- compare.mechanisms(comparison)
    filename = sprintf("out/%s%d-m%dn%d-results.Rdata",
                       filename.prefix,
                       ifelse(include.3way, 3, 2), nhospitals, nsize)
    cat(sprintf("\nSaving to filename %s", filename))
    save(results, file=filename)
    return(results)
  }
  exp.results <- list()
  if(1 %in% run.profiles)
    exp.results = run.comparison(exp.results, nhospitals, nhospitals-1, T)
  if(2 %in% run.profiles)
    exp.results = run.comparison(exp.results, nhospitals, nhospitals-1, F)
  if(3 %in% run.profiles)
    exp.results = run.comparison(exp.results, 1, 0, T)
  if(4 %in% run.profiles)
    exp.results = run.comparison(exp.results, 1, 0, F)
}

run.all.simple.experiments <- function(nsamples=1000) {
  Simple1 <- simple.experiments(1, nsamples=nsamples, max.hospitalSize=150)
  Simple2 <- simple.experiments(2, nsamples=nsamples, max.hospitalSize=150)
  results = list(simple1=Simple1,
                 simple2=Simple2)
  save(results, file="out/SimpleExperiments.Rdata")
  print("Done with Simple Experiments")
}

simple.experiments <- function(experiment.no, nsamples=100, max.hospitalSize=140) {
  ## Last suite of simpler experiments
  results = list(UPRA=matrix(0, nrow=0, ncol=3), 
                 NonUPRA=matrix(0, nrow=0, ncol=3))
  print(sprintf("Experiment id=%d Samples=%d Max.Hospitalsize=%d",
                experiment.no, nsamples, max.hospitalSize))
  add.result <- function(pra.value, results, update) {
    key = ifelse(pra.value, "UPRA", "NonUPRA")
    CHECK_EQ(length(update), ncol(results$UPRA))
    results[[key]] <- rbind(results[[key]], update)
    colnames(results[[key]]) <- c("nsize", "total", "matched")
    rownames(results[[key]]) <- NULL
    return(results)
  }
  
  if(experiment.no == 1) {
    nsize.list = as.integer(seq(20, max.hospitalSize, by=25))
    pra.list = c(T, F)
  
    add.result <- function(pra.value, results, update) {
      key = ifelse(pra.value, "UPRA", "NonUPRA")
      CHECK_EQ(length(update), ncol(results$UPRA))
      results[[key]] <- rbind(results[[key]], update)
      colnames(results[[key]]) <- c("nsize", "total", "matched")
      rownames(results[[key]]) <- NULL
      return(results)
    }
    
    pb = txtProgressBar(style=3)
    
    for(pra in pra.list) {
      for(i in 1:nsamples) {
        nsize = sample(nsize.list, size=1, replace=T)
        # print(sprintf("pra=%s n=%d", pra, nsize))
        rke = rrke(nsize, uniform.pra=pra)
        xRKE = rke.extended.Rsubgraph(rke)
        Rpair.ids = subset(xRKE$pairs, pair.type=="R")$pair.id
        m = max.matching(xRKE, include.3way=T, promote.pair.ids=Rpair.ids)
        Rpair.ids.matched = subset(m$match, pair.type=="R")$pair.id
        results <- add.result(pra, results, c(nsize, length(Rpair.ids), length(Rpair.ids.matched)))
        setTxtProgressBar(pb, value = nsamples * (1-pra) + i / ( nsamples))
      }
    }
    
    return(results)
  } else if(experiment.no==2) {
    nHospitals = 3  # total no. of hospitals
    nsize.list = as.integer(seq(20, max.hospitalSize, by=25))
    pb = txtProgressBar(style=3)
    for(i in 1:nsamples) {
      nsize = sample(nsize.list, size=1, replace=T)
      pool = rrke.pool(m=nHospitals, n=nsize, uniform.pra=T)
      all.c = paste(rep("c", nHospitals), collapse="")
      kpd = kpd.create(pool, strategy.str=all.c, include.3way=F)
      rke.reported = kpd$reported.pool$rke.all
      nAB = subset(rke.reported$pairs, desc=="A-B")$pair.id
      nBA = subset(rke.reported$pairs, desc=="B-A")$pair.id
      
      m = max.matching(rke.reported, promote.pair.ids=c(nAB, nBA))
      nAB.matched = subset(m$match, desc=="A-B")$pair.id
      nBA.matched = subset(m$match, desc=="B-A")$pair.id
      total.R = length(nAB) + length(nBA)
      total.matched = length(nAB.matched) + length(nBA.matched)
      ## add result
      results <- add.result(T, results, c(nsize, total.R, total.matched))
      setTxtProgressBar(pb, value = i / ( nsamples))
    }
    return(results)
  } else {
    stop(sprintf("Simple Experiment id=%d not valid", experiment.no))
  }
  
}

run.sweetSpot.experiments <- function(nsamples) {
  table.welfare.incentives(mechanisms=c("xCM"), nhospitals=3, nsize=160, 
                           include.3way=F, nsamples=nsamples,
                           filename.prefix="SweetSpot", full.suite=F)
  
  table.welfare.incentives(mechanisms=c("xCM"), nhospitals=2, nsize=80, 
                           include.3way=T, nsamples=nsamples,
                           filename.prefix="SweetSpot", full.suite=F)
}

