##  Tables of the paper. 
## Code to save into latex table.
library(xtable)

# Add necessary libs.
source("../r-toolkit/checks.R")
source("terminology.R")
source("rke.R")
source("matching.R")
source("mechanisms.R")
## hand-made bootstrap (USDA organic)
bootstrap = function(x) sd(replicate(1000, { y = sample(x, replace=T); mean(y)}))

## Combines the result matrix with the SE errors.
add.se = function(M, SE) {
    D = matrix(NA, nrow=nrow(M), ncol=ncol(M))
    colnames(D) = colnames(M)
    for(i in 1:nrow(M)) {
        D[i,]= sapply(1:ncol(M), function(j) sprintf("%.2f (%.2f)", M[i,j], SE[i,j]))
    }
    return(D)
}

compare.mechanisms.kpd <- function(mechanisms,
                                   m, include.3way,
                                   kpd.baseline,
                                   kpd.deviation) {
  # Compares two mechanisms.
  # m = #hospitals
  # Result looks like this:
  #   baseline=> { strategy="...", 
  #                rCM=LIST(utility=MATRIX(m x Ntrials), matching.info),
  #                ...}
  #   deviation = > {...}
  CHECK_MEMBER(mechanisms, c("rCM", "xCM", "Bonus"))
  result = list(baseline=list(strategy="na"),
                deviation=list(strategy="na"))
  # Initialize
  for (mech in mechanisms) {
    for (j in names(result))
      result[[j]][[mech]] <- list(utility=rep(0, m),
                                  matching.info=empty.match.result(empty.rke())$information)
  }
  
  for (mech in mechanisms) {
    m.baseline = Run.Mechanism(kpd=kpd.baseline, mech=mech, include.3way=include.3way)$total.matching
    m.deviation = Run.Mechanism(kpd=kpd.deviation, mech=mech, include.3way=include.3way)$total.matching
    
    result$baseline[[mech]]$utility <- get.matching.hospital.utilities(m.baseline, m)
    result$baseline[[mech]]$matching.info = m.baseline$information
    result$baseline[[mech]]$matching = m.baseline
        
    result$deviation[[mech]]$utility <- get.matching.hospital.utilities(m.deviation, m)
    result$deviation[[mech]]$matching.info = m.deviation$information
    result$deviation[[mech]]$matching = m.deviation
  }
  return(result)                       
}



compare.mechanisms.strategy <- function(mechanisms,
                               baseline.strategy,
                               deviation.strategy,
                               m, n, include.3way, uniform.pra,
                               ntrials=10) {
  # Compares two mechanisms.
  # Result looks like this:
  #   baseline=> { strategy="...", 
  #                rCM=LIST(utility=MATRIX(m x Ntrials), matching.info),
  #                ...}
  #   deviation = > {...}
  CHECK_MEMBER(mechanisms, c("rCM", "xCM", "Bonus"))
  result = list(baseline=list(strategy=baseline.strategy),
                deviation=list(strategy=deviation.strategy))
  # Initialize
  for (mech in mechanisms) {
    for (j in names(result))
      result[[j]][[mech]] <- list(utility=matrix(0, nrow=m, ncol=ntrials),
                                  matching.info=empty.match.result(empty.rke())$information)
  }
  
  pb = txtProgressBar(style=3)
  for(i in 1:ntrials) {
    rke.pool = rrke.pool(m=m, n=n, uniform.pra=uniform.pra)
    kpd.baseline = kpd.create(rke.pool=rke.pool, strategy.str=baseline.strategy)
    kpd.deviation = kpd.create(rke.pool=rke.pool, strategy.str=deviation.strategy)
    
    single.result <- compare.mechanisms.kpd(mechanisms,m=m, include.3way=include.3way,
                                            kpd.baseline=kpd.baseline,
                                            kpd.deviation=kpd.deviation)
    for (mech in mechanisms) {
      result$baseline[[mech]]$utility[, i] <- single.result$baseline[[mech]]$utility
      result$baseline[[mech]]$matching.info = result$baseline[[mech]]$matching.info + 
        single.result$baseline[[mech]]$matching.info
      
      
      result$deviation[[mech]]$utility[, i] <- single.result$deviation[[mech]]$utility
      result$deviation[[mech]]$matching.info = result$deviation[[mech]]$matching.info +
        single.result$deviation[[mech]]$matching.info
    }
    setTxtProgressBar(pb, value=i/ntrials)
  }
  return(result)                       
}


SingleDeviation.experiments <- function(deviate.to="r", mech,
                                        m=4, n=20, ntrials=1,
                                        verbose=F) {
  # These experiments check the utility gains 
  # of one hospital going from "truthful" to "deviate.to" strategy.
  # We take "ntrials" samples and run on "m" hospitals with "n" pairs each.
  #
  pb = txtProgressBar(style=3)
  performance.deviator <- c()
  
  for(trial in 1:ntrials) {
    # 1. Sample RKE pool
    rp <- rrke.pool(m=m, n=n, uniform.pra=T)
    
    # Which hospital is going to deviate?
    deviate.hid <- sample(1:m, size=1)
    deviate.factor <- 0
    # If deviation is "R" then make sure it makes sense, i.e.
    # pick the hospital for which we expect better results.
    # 
    # 2. Compute A-B/B-A pairs/hospital
    # col1-2 = AB/BA of hospital,   Cols3-4=AB/BA rest of hospitals
    # col5 = R-attack successful
    # 
    if (deviate.to=="r") {
      ABpairs <- matrix(0, nrow=m, ncol=5)
      # Compute the A-B pairs for each hospitals
      for (i in 1:m) {
        ABpairs[i, 1] <- nrow(subset(rp$rke.all$pairs, desc=="A-B" & hospital==i))
        ABpairs[i, 2] <- nrow(subset(rp$rke.all$pairs, desc=="B-A" & hospital==i))
        ABpairs[i, 3] <- nrow(subset(rp$rke.all$pairs, desc=="A-B" & hospital!=i))
        ABpairs[i, 4] <- nrow(subset(rp$rke.all$pairs, desc=="B-A" & hospital!=i))
        longSide = which.max(c(ABpairs[i, c(1,2)]))  # 1=A-B, 2=B-A
        shortSide = setdiff(c(1,2), c(longSide))
        otherOppSide = 5 - longSide  # e.g. 1(A-B) -> 4 (B-A), 2->3
        otherSameSide = longSide + 2  # e.g 1=A-B -> 3 (A-B)
        
        expected.matches <- ABpairs[i, longSide] * ABpairs[i, otherOppSide] /
          (ABpairs[i, longSide] + ABpairs[i, otherSameSide])
        surplus = abs(ABpairs[i,1] - ABpairs[i,2])
        # Expected gain = surplus - expected.matches
        # if E[matches] > surplus this is bad. 
        # E[matches] < surplus  is unlikely to happen (hospital sizes)  
        # Ideally we want to match our surplus exactly
        # So, if the below is negative, then this is bad news.
        #   and the closer to the 0 the better.
        ABpairs[i, 5] <- surplus - expected.matches
      }
      
      colnames(ABpairs) <- c("A-B(H)", "B-A(H)", "A-B(rest)", "B-A(rest)", "Gain-R-Dev")
      # 3. Choose which hospital will be deviating
      #    The one where the surplus is best matched.
      deviate.hid = which.max(ABpairs[, 5])
      deviate.factor = max(ABpairs[, 5])
      
      if(verbose) {
        print("A-B pairs are")
        print(ABpairs)
       }
    }
    # Create strategy
    strategy = rep("t", m)
    strategy[deviate.hid] <- deviate.to
    if(deviate.to == "r") {
      if (deviate.factor < 0.2) 
        strategy[deviate.hid] <- "t"
    }
    # Define strategies.
    # baseline = "ttttt..."
    # deviation = "tt...x...ttt"  where x=deviate.to
    base.strategy <- paste(rep("t", m), collapse="")
    dev.strategy <- paste(strategy, collapse="")
    
    # Create the baseline/deviation KPDs and then compare mechanisms.
    kpd.dev <- kpd.create(rke.pool=rp, strategy.str=dev.strategy)
    kpd.base <- kpd.create(rke.pool=rp, strategy.str=base.strategy)
    # Output values.
    x = compare.mechanisms.kpd(mechanisms=c(mech), m=m, include.3way=F,
                               kpd.baseline=kpd.base, kpd.deviation=kpd.dev)
    x$deviate.hid = deviate.hid
    setTxtProgressBar(pb, value=trial/ntrials)
    # Save how well the deviator improved.
    performance.deviator <- c(performance.deviator, 
                              x$deviation[[mech]]$utility[deviate.hid] -  x$baseline[[mech]]$utility[deviate.hid])
    
    if(verbose) {
      print(sprintf("Deviating hospital=%d", deviate.hid))
      print(sprintf("Strategies base=%s dev=%s", base.strategy, dev.strategy))
      
      print(sprintf("%s before:%d -> %d", mech,
                    x$baseline[[mech]]$utility[deviate.hid],
                    x$deviation[[mech]]$utility[deviate.hid]))
      print(sprintf("mu=%.3f  se=%.3f", mean(performance.deviator),
                    bootstrap.mean(performance.deviator)))
    }
   
    save(performance.deviator, file="out/Rdev-experiment.Rdata")
  }
}

