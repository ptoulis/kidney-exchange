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

summarize.coompare.output <- function(compare.mech.out) {
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

table1.theoretical.violations <- function(nsamples=100) {
  # This experiment explores two things
  # 1. The μ(n) formula = expected #matches in Gn
  # 2. The PM assumption, by checking on the matches of regular matching.
  #
  # Table1 is of the following form:
  #  n (size)  #matches
  # 
  # Table 2 is of the form
  # n(size)  (OR OO OS....)  = MATCHING INFO
  
  all.sizes <- round(seq(5, 200, by=5))
  sampled.sizes <- sample(all.sizes, size=nsamples, replace=T)
  print("Sampled sizes breakdown")
  print(table(sampled.sizes))
  # stores #matches for each RKE
  nmatches <- c()
  # stores #R pairs in each RKE (this is useful for the theoretical #matches)
  nRpairs <- c()
  # stores information about the matches e.g. OR, OO, OS, ...
  match.information = matrix(NA, nrow=0, ncol=length(empty.match.result(empty.rke())$information))
  # format of the output object.
  get.result.object <- function() {
    CHECK_EQ(length(nRpairs), length(nmatches))
    CHECK_EQ(length(nRpairs), nrow(match.information))
    return(data.frame(nsize=head(sampled.sizes, length(nRpairs)), 
               Rsize=nRpairs,
               matched=nmatches, 
               match=match.information))
  }
  pb = txtProgressBar(style=3)
  CHECK_EQ(length(sampled.sizes), nsamples)
  
  # save at the last iteration too.
  save.checkpoints <- c(seq(1, nsamples-1, length.out=20), nsamples)
  
  for(i in 1:nsamples) {
    n = sampled.sizes[i]
    rke = rrke(n)
    m = max.matching(rke, include.3way=F, regular.matching=T)
    nmatches <- c(nmatches, get.matching.utility(m))
    nRpairs <- c(nRpairs, nrow(subset(rke$pairs, pair.type=="R")))
    match.information <- rbind(match.information, m$information)
    setTxtProgressBar(pb, value=i/nsamples)
    
    if(i %in% save.checkpoints) {
      table1 = get.result.object()
      save(table1,file="out/table1.Rdata")
    }
  }
  print("Simulation complete. File saved in out/table1.Rdata")
}

table2.welfare.incentives.2way <- function(nsamples=100) {
  # Table of 2way exchanges to compare Welfare and Incentives.
  #
  results = list()
  nhospitals = 6
  get.strategy.profile <- function(no.truthful) {
    CHECK_TRUE(no.truthful >= 0 & no.truthful <= nhospitals, msg="#truthful should be correct")
    no.deviating = nhospitals - no.truthful
    return(paste(c(rep("t", no.truthful), rep("c", no.deviating)), collapse=""))
  }
  
  run.experiment <- function(base.Nt, dev.Nt, pra) {
    # Runs a single experiment.
    # 
    # Args:
    #  base.Nt = no. of truthful agents in baseline strategy profile
    #  dev.Nt = no. of truthful strategies in dev strategy profile
    #  pra = T or F, whether we want uniform PRA or non-uniform PRA.
    baseline.strategy = get.strategy.profile(base.Nt)
    deviation.strategy = get.strategy.profile(dev.Nt)
    print(sprintf("Comparing profiles  %s vs. %s", baseline.strategy, deviation.strategy))
    comparison = create.comparison(mechanisms=c("rCM", "xCM", "Bonus"),
                                   nHospitals=nhospitals, nSize=20,
                                   uniform.pra=pra, 
                                   include.3way=F,
                                   baseline.strategy=baseline.strategy,
                                   deviation.strategy=deviation.strategy,
                                   nsamples=nsamples)
    profile.name = sprintf("prof%d%d", base.Nt, dev.Nt)
    result.name = sprintf("%s-%s", profile.name, ifelse(pra, "UPRA", "NonUPRA"))
    results[[result.name]] <- compare.mechanisms(comparison)
    save(results, file="out/table2-results.Rdata")
  }
  
  run.experiments(nhospitals, nhospitals-1, T)
  run.experiments(nhospitals, nhospitals-1, F)
  run.experiments(1, 0, T)
  run.experiments(1, 0, F)  
}





