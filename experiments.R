##  Tables of the paper. 
## Code to save into latex table.
library(xtable)

rm(list=ls())
# Add necessary libs.
source("../r-toolkit/checks.R")
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
  # For example,
  #   result$baseline$rCM$utility = (m x samples) utility matrix
  #   result$baseline$rCM$internal.matchingInfo[[3]] = MATCHING_INFO for hospital 3
  #   result$baseline$rCM$mech.matchingInfo = MATCHING INFO for only the mechanism
  all.strategies <- c("baseline", "deviation")
  result = list(baseline=list(strategy=comparison$baseline.strategy),
                deviation=list(strategy=comparison$deviation.strategy))
  
  empty.matchingInfo <- empty.match.result(empty.rke())$information
  # 0. Initialize the result object.
  for (mech in comparison$mechanisms) {
    for (s in all.strategies) {
      # Setting the result object for (strategy, mechanism)
      result[[s]][[mech]] <- list(utility=matrix(0, nrow=comparison$m, ncol=comparison$nsamples),
                                  mech.matchingInfo=empty.matchingInfo,
                                  internal.matchingInfo=list(),
                                  total.matchingInfo=empty.matchingInfo)
      
      for (hid in 1:comparison$m) {
        # Initialize the matching info
        result[[s]][[mech]]$internal.matchingInfo[[hid]] = empty.matchingInfo
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
        }  # for all hospitals, matching info
      } # all strategies (baseline, deviation)
      
    } # for every mechanism
    setTxtProgressBar(pb, value=i / comparison$nsamples)
  }  # all trials
  
  ## Return the final result
  return(result)                       
}
