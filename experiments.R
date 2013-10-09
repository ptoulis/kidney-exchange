##  Tables of the paper. 
## Code to save into latex table.
library(xtable)
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

## Dumps the table as a .tex file.
table.to.tex = function(results, filename) {
    fileConn = file(filename,open="w")
    lines =c()
    attachLine <- function(s) {
        lines <<- c(lines, s)
    }
    attachLine("\\begin{center}")
    
    columns = colnames(results)
    ncols = length(columns)
    ##    alignment
    attachLine(sprintf("%s%s}","\\begin{tabular}{", 
                       paste(rep("c", ncols),  collapse=" | ")))
    tHeader = c(columns)
    attachLine(sprintf("%s\\\\ \\hline", paste(tHeader, collapse=" & ") ))
    for(i in  1:nrow(results)) {
        values = results[i,]
        values = sapply(values, function(x) ifelse(is.numeric(x), sprintf("%.2f", x), x))
        attachLine(sprintf("%s \\\\",paste(values, collapse=" & ")))
    }
    attachLine("\\end{tabular}")
    attachLine("\\end{center}")
    writeLines(lines,con=fileConn)
    close(fileConn)
}

tables.all.12 <- function() {
  # Table 1. Square-root law.
  loginfo("Generating Table 1: mu(n) theorem and assumptions")
  D = table.matchings(sizes = c(50, 100, 200, 300), m=2, trials=1000);
  table.to.tex(D, filename="out/tex/table1-assumptions.tex");
  
  # Table 2. tab
  loginfo("Generating Table 2: Violations.")
  D = table.violations(sizes = c(50, 100, 200, 300), trials=1000);
  table.to.tex(D, filename="out/tex/table2-violations.tex");
}

# Run all tests at once. 
tables.all = function(sizes.array=c(20), ntrials=200) {
    
  # Tables 3,4,5
  loginfo("Generating Scenarios table for rCM")
  table.mechs("rCM", m=3, sizes=sizes.array,
              include.3way=T, trials=ntrials);
  loginfo("Generating Scenarios table for xCM")
  table.mechs("xCM", m=3, sizes=sizes.array, 
              include.3way=T, trials=ntrials);
  loginfo("Generating Scenarios table for Bonus")
  table.mechs("Bonus", m=4, sizes=sizes.array, 
              include.3way=T, trials=ntrials);
  table.mechs.to.graph();
  
  # Table 6, 7
  loginfo("Generating efficiency table with uniform PRA.")
  table.efficiency(m=4, sizes=sizes.array,
                   include.3way=T, uniform.pra=T, trials=ntrials,
                   filedesc="efficiency");
  loginfo("Generating efficiency table with non-uniform PRA.")
  table.efficiency(m=4, sizes=sizes.array, 
                   include.3way=T, uniform.pra=F, trials=ntrials,
                   filedesc="efficiency");
  table.efficiency.to.graph()
  return()
  # Table 8
  loginfo("Generating efficiency table with varying m (#of hospitals)")
  table.efficiency.many.hospitals(many.m=c(4,6,8,10,14,18,20), n=25, trials=200)
  table.efficiency.many.to.graph()
}

###      FUNCTIONS TO CREATE THE TABLES
##  Table 1.   Theorem 1 (expected value of maximum matchings.)
##    n | selfish | together | m(n)   | surplus
##   Synopsis.    D = table.matchings(sizes = c(50,100,200, 300), m=2, trials=100)
##               table.to.tex(D, filename="out/assumptions.tex")
table.matchings = function(sizes=c(50), 
                           m= 3, 
                           trials=10) {
  colnames = c("size ($n$)", "selfish", "together", "together (theoretical)", "surplus/$\\sqrt{n}$")
  ncols = length(colnames) ## size-experimental matches - theoretical matches  - (same for non-PRA)
  pb = txtProgressBar(style=3)
  ## result matrix
  M = matrix(NA, nrow=length(sizes), ncol=ncols)
  ## standard error matrix
  SE =  matrix(NA, nrow=length(sizes), ncol=ncols)
  colnames(M) = colnames
  ##  Total number of runs
  nruns = length(sizes) * trials
  for(i in 1:length(sizes)) {
    n = sizes[i]  # hospital size
    ##  x =  2 x trials
    x = sapply(1:trials, function(j) {
      val = (i-1) * trials + j
      setTxtProgressBar(pb, val/nruns)
      
      # 0. Sample rke's 
      rke.pool = rrke.pool(m=m, n=n, uniform.pra=T)
      rke.list = rke.pool$rke.list
      CHECK_EQ(length(rke.list), m)
      # 1.   Selfish (each one matches in isolation) 
      selfish.matches = c()
      for(hid in 1:length(rke.list)) {
        selfish.matches[hid] = get.matching.utility(max.matching(rke.list[[hid]]))
      }
      selfish.matches = sum(selfish.matches)
      
      # 2. Pool and match together
      rke.all = rke.pool$rke.all
      m.all = max.matching(rke.all)
      together.matches = get.matching.utility(m.all)
      
      # 3. Theoretical matches
      together.theoretical = mu.thm(n=n, m=m)
      
      # 4. Surplus
      surplus = (together.matches - selfish.matches)/sqrt(n)
      
      return(c(selfish.matches, together.matches, together.theoretical, surplus))
    });
    
    M[i,] = c(n,  t(apply(x, 1, mean)) )
    SE[i,] = c(0, t(apply(x, 1, bootstrap) ))
  }   
  # text = 
  D = add.se(M, SE)
  print("")
  return(D)
}


##  Table 2. Violations
##  Sample usage:
##  ------------------------------------------------------
##  D = table.violations(c(50, 100..),  sims=100)
##  table.to.tex(D, filename="out/violations.tex")
table.violations = function(sizes=c(50), trials=10) {
     x = E2(sizes=sizes, sims=trials)
     ##   x[[size]]
     ##  Results matrix.
     M = matrix(NA, nrow=length(sizes), ncol=7)
     SE = matrix(NA, nrow=length(sizes), ncol=7)
     for(i in 1:length(sizes)) {
         result = x[[i]]
         M[i,] = apply(result, 2, mean)
         SE[i,] = apply(result, 2, bootstrap)
     }
     D = add.se(M, SE)
     return(D)
}
##   
##     Experiments  2 :   Violations of Regularity and PM assumptions.
##
###   results[[i]] = sims x 6

##    i = hospital size (=20 * i)
##    sims = no. of simulations
##   The data columns are: Regular,  R,  O-O, A-A, B-B, AB-AB
##     where each one counts the no. of violations for each case (regularity, R, and S subgraphs)
##   1.  To get the mean regularity violations at size=60
##         apply(res[[3]], 2, mean)[,1]    (b/c column 1 has regularity violations)
E2 <- function(sizes=c(50), sims=100) {
  col.names= colnames(e2.once(20))
  
  results = list()
  pb = txtProgressBar(style=3)
  total.N = length(sizes) * sims
  val = 0
  for(i in 1:length(sizes)) {
    size = sizes[i]
    results[[i]] = matrix(0, nrow=0, ncol=length(col.names))
    colnames(results[[i]]) = col.names
    for(j in 1:sims) {
      ## Calling e2.once()
      res = as.vector(e2.once(size))
      results[[i]] = rbind(results[[i]], c(res))
      val = val+1
      
      setTxtProgressBar(pb, value=val/total.N)
    }
  }
  return(results)
}
## Violations (for one specific size)
e2.once <- function(size) {
  #ns = 20 * 1:size
  #A = matrix(0, a)
  
  ## Sample the RKE object 
  rke = rrke(size)
  
  ## Get OR, OS edges.
  OR = filter.edges.by.type(rke,"O", "R")
  OS = filter.edges.by.type(rke, "O","S")
  to.remove =union(OR,OS)
  
  ###
  Rs = c()
  ## Return # of A-B or B-A pairs
  getRtypeL <- function(arg_rke, Rtype) {
    if(Rtype==1)
      return(length(filter.pairs.by.donor.patient(arg_rke, "A","B")))
    return(length(filter.pairs.by.donor.patient(arg_rke, "B","A")))
  }
  Rs[1] = getRtypeL(rke, 1)
  Rs[2] = getRtypeL(rke, 2)
  r.short.side = which.min(Rs)
  ## Do the same or S
  
  m = max.matching(rke, remove.edges=to.remove)
  
  rke2.remainder = remove.pairs(rke, m$matching$matched.ids)
  types.left = get.pairs.attribute(rke2.remainder, "type")
  
  ## 1. Regularity violations = count remaining OD's
  RegViols = sum(types.left=="O")
  if(RegViols>3)
    save(rke, file="debug/e2.Rdata")
  ## 2. R-violations = count remaining size of short R side
  RViols = getRtypeL(rke2.remainder, r.short.side)
  ## 3. S-violations = count >2 for each X-X
  SViols = sapply(Blood.Types, function(abo) 
    max(0,length(filter.pairs.by.donor.patient(rke2.remainder, abo,abo))-1))
  ## Done checking violations
  M  = matrix(c(size, RegViols, RViols, SViols), nrow=1, ncol=7)
  colnames(M) = c("n","Reg.", "R", "OO", "AA","BB","ABAB")
  return(M)
}
## Compute the relative gain for _hid.interest_ 
# mechanism _mech_ under _setup_, for two different strategies
# Returns:  2 x trials   matrix of utiliies.
relative.gain = function(kpd1, kpd2, mech, hid, include.3way) {
  # Runs both KPD markets for the same mechanism
  # Then returns the utility for the particular hospital (hid)
  #
  # Returns:
  #   A 1x2 vector of two utilities
  U1 = Run.Mechanism(kpd=kpd1, mech=mech, include.3way=include.3way)
  U2 = Run.Mechanism(kpd=kpd2, mech=mech, include.3way=include.3way)
  c(U1[hid], U2[hid])
}

relative.gain.scenario = function(scenario, mech, m, n, include.3way, trials,
                                  pb=txtProgressBar(max=trials, style=3),
                                  pb.start=0, verbose=F) {
  # Runs the scenario for the particular mechanism
  #
  # Args:
  #   scenario = one of {A, B, C, D} as a character
  #   mech = mechanism (rCM, xCM, Bonus)
  #   m = #hospitals
  #   n = size of hospitals
  #   trials= #samples to take.
  #   pb = progress bar
  #   pb.start = where to start the progress bar from
  #
  # Returns:
  #   A 2 x trials matrix that has utilities.
  #   Row 1 has the utility of the strategic H1
  #   Row 2 has the utility for the baseline "ttt..." scenario
  #
  # What H-1 is doing per mechanism, per scenario
  h1.str.list = list("rCM" = list(A="c", B="c", C="c", D="c"),
                     "xCM" = list(A="c", B="r", C="c", D="r"),
                     "Bonus"=list(A="c", B="r", C="c", D="r") )
  # What other hospitals are doing per mechanism, per scenario
  others.list = list("rCM" = list(A="t", B="c", C="t", D="c"),
                     "xCM" = list(A="t", B="t", C="t", D="t"),
                     "Bonus"=list(A="t", B="t", C="t", D="t"))
  
  uniform.pra.list = list(A=T, B=T, C=F, D=F)
  
  ## define the parameters based on mechanism + scenario
  h1.strategy = h1.str.list[[mech]][[scenario]]
  others = rep(others.list[[mech]][[scenario]],  m-1)
  is.uniform.pra = uniform.pra.list[[scenario]]  
  # output matrix
  A = matrix(NA, nrow=2, ncol=trials)
  ## e.g. Scenario A ->  str1 = "ctt" str2 = "ttt"
  str1 = paste(c(h1.strategy, others), collapse="")
  ## str2 =  baseline strategy profile.
  str2 = paste(c("t", others), collapse="")
  if(verbose) {
    cat("\n")
    loginfo(sprintf("Scenario %s, Mechanism %s TEST, size=%d, 3-chains=%s: profile %s vs. baseline %s, uniform.pra=%s",
                  scenario, mech, n, include.3way, str1, str2, is.uniform.pra))
  }
  rownames(A) <- c(sprintf("%s(test)", str1), sprintf("%s(baseline)", str2))
  for(i in 1:trials) {
    rke.pool = rrke.pool(m=m, n=n, uniform.pra=is.uniform.pra)
    # Create the Kidney-Paired Donation market
    kpd1 = kpd.create(rke.pool, str1)
    kpd2 = kpd.create(rke.pool, str2)
    
    utils = relative.gain(kpd1, kpd2, mech=mech, hid=1, include.3way=include.3way)
    A[1,i] = utils[1] # has the utility for KPD1
    A[2,i] = utils[2]  # has the utility for KPD2
    setTxtProgressBar(pb, value=pb.start + i)
  }
  if(verbose) {
    cat("\nSimulation over")
  }
  return(A)
}

table.mechs = function(mech, m=3, sizes=c(20), include.3way, trials=10, verbose=F) {
  # Creates Tables 3,4,5. Runs all scenarios using relative.gain.scenario
  # and all different sizes.
  # Sample usage:
  #   table.mechs("rCM", m=3, sizes=c(20, 40), trials=100)
  #   -> TABLE 3 -> saves a file as "out/mech-rCM-m3-trials100.Rdata"
  #   table.mechs("xCM", m=3, sizes=...., trials=..) 
  #   -> TABLE 4  (etc)
  # Args:
  #   mech = mechanism (character) i.e. "rCM", "xCM" etc
  #   m = #hospitals
  #   sizes = vector of hospital sizes
  #   trials = #samples for each size
  # Returns:
  #   results[[scenario.id]][[size]] = Matrix (relative.gain.scenario output)
  results = list()  
  ## each mechanism has 4 scenarios
  scenarios = c("A", "B", "C", "D")
  print(sprintf("Running scenarios for mech=%s m=%d #sizes=%d", mech, m, length(sizes)))
  # total number of trials
  N = length(scenarios) * length(sizes) * trials
  pb = txtProgressBar(style=3,min=0, max=N)
  # output filename
  filename = sprintf("out/mech-%s-m%d-trials%d.Rdata", mech, m, trials)
  for(sce.i in 1:length(scenarios))   {
    scen = scenarios[sce.i]  # load the scenario
    results[[scen]] = list()
    for(size.j in 1:length(sizes)) {
      n = sizes[size.j]
      pb.start = length(sizes) * trials * (sce.i-1) + trials * (size.j-1)
      results[[scen]][[sprintf("%d",n)]] = relative.gain.scenario(scenario=scen,
                                                                  mech=mech,
                                                                  m=m, n=n,
                                                                  include.3way=include.3way,
                                                                  trials=trials,
                                                                  pb=pb, pb.start=pb.start,
                                                                  verbose=verbose)
      save(results, file=filename)
    }
  }
  return(results)
}

table.mechs.to.graph = function() {
  # Traverses the "out/" folder for "mech-*" files
  # loads and then saves results in a .png graph
  files = list.files(path="out/", pattern="mech", full.names=T)
  for(filename in files) {
    x = filename
    m = regexec(text=x, pattern="mech-(.*?)-")
    mech = regmatches(x, m)[[1]][2]
    load(filename)
    # assume:  results.
    png.file = sprintf("out/png/tables345-mech-%s.png", mech)
    loginfo(sprintf("Mechanism %s, Saving to %s ", mech, png.file))
    png(file= png.file)
    par(mfrow=c(2,2) )
    for(i in names(results) ) {
      obj=  results[[i]]
      for(j in names(obj)) {
        A = obj[[j]]  ## 2 x trials 
        obj[[j]] = A[1,] / A[2,]
      }
      boxplot(obj, ylim=c(0.5, 2), main=sprintf("mech=%s scenario= %s", mech, i ), ylab="ratio (str/base)", xlab="size")
    }
    dev.off();
  }
}

table.efficiency = function(m=4, sizes=c(20), uniform.pra, include.3way,
                            trials=10, filedesc="efficiency") {
  # Creates Table 6 and 7
  # Sample usage
  # -------------------------------------
  # > table.efficiency(m=4, sizes=c(20, 40, 60, 80), uniform.pra=T, trials=100, 
  #                 filedesc="efficiency")
  # saves ONE file in "out/efficiency-pra-TRUE-m4-trials100"
  # > table.efficiency.to.graph() 
  results = list()
  # Compare rCM at canonical deviation
  # xCM at truthful
  # Bonus at r-strategy
  # all.str("t") = "ttt..." length = #hospitals (m)
  all.str = function(ch)  paste(rep(ch, m), collapse="")
  
  # config = configuration of the experiment
  config = list("xCM_t" = list(mech="xCM", str=all.str("t")),
                "rCM_c" = list(mech="rCM", str=all.str("c")),
                "Bonus_c" = list(mech="Bonus", str= all.str("c"))) # this seems to hurt Bonus
                # "Bonus_r" = list(mech="Bonus", str=all.str("r")))
  # N = total #iterations
  N = length(sizes) * trials * length(names(config))
  # progress bar (i really like them)
  pb = txtProgressBar(style=3, min=0, max=N)
  count = 0
  for(i in 1:length(sizes)) {
    # i = index of size
    # n = size of hospital
    n = sizes[i]
    n.str = sprintf("%d", n)
    # Initialize the mechanisms
    for(mech.name in names(config) ) {
      results[[mech.name]] = list()
      results[[mech.name]][[n.str]] = c()
    }
    cat("\n")
    loginfo(sprintf("Efficiency m=%d, uniform=%s, size=%d, 3-chain=%s", m, uniform.pra, n, include.3way));
    for(j in 1:trials) {
      # j = trial
      # Sample rke.list/rke.all
      rke.pool = rrke.pool(m=m, n=n, uniform.pra=uniform.pra)
      # Iterate over mechanisms to figure out efficiency
      for(mech.name in names(config)) {
        # (i, j, mech.name) = (size, trial, mechanism)
        strategy = config[[mech.name]][["str"]]
        mech = config[[mech.name]][["mech"]]
        kpd = kpd.create(rke.pool, strategy)
        welfare = sum(Run.Mechanism(kpd=kpd, mech=mech, include.3way=include.3way))
        # Append welfare to the mechanism vector.
        results[[mech.name]][[n.str]] = c( results[[mech.name]][[n.str]], welfare )
        count = count + 1
        setTxtProgressBar(pb, value=count)
      } # iterate over all mechanisms.
      
      save(results, file=sprintf("out/%s-pra-%s-m%d-trials%d.Rdata", 
                                 filedesc, uniform.pra, m, trials))
    }  # trials
  } # iterate over all sizes.
  
  return(results)
}

# Takes all files in "experiments" that match "efficiency-" and 
# creates a 2x2 boxplot, where block = mechanism, x-axis = hospital size, 
# y-size = efficiency ratio
table.efficiency.to.graph = function() {
  files = list.files(path="out/", pattern="efficiency-", full.names=T)
  for(filename in files) {
    x = filename
    m = regexec(text=x, pattern="pra-(.*?)\\.")
    uniform.pra = regmatches(x, m)[[1]][2]
    load(filename)
    # assume:  results.
    png.file = sprintf("out/png/efficiency-uniform-PRA-%s.png", uniform.pra)
    cat(sprintf("\nData filename: %s Mechanisms = %s", filename, paste(names(results), collapse=", ")))
    cat(sprintf("\nImage filename %s \n", png.file))
    png(file= png.file)
    par(mfrow=c(2,2) )
    for(mech in names(results) ) {
      obj = results[[mech]]
      for(size in names(obj)) {
        A = obj[[size]]  ## 2 x trials 
        A2 = results[["rCM_c"]][[size]]
        obj[[size]] = A / A2
        cat(sprintf("\nRelative efficiency for mech=%s is %.3f\n", mech, mean(obj[[size]])))
        cat("Summary:\n")
        print(summary(obj[[size]]))
      }
      boxplot(obj, ylim=c(0.6, 1.1), main=sprintf("mech=%s", mech),
              ylab="ratio (str/base)", xlab="size")
    }
    dev.off();
  }
}

## Efficiency for many hospitals
# Table 8: Sample usage
# ------------------------
# table.efficiency.many.hospitals(many.m = c(4, 6, 8, 10, 14, 16), 
#                                 n=25, trials=100, filedesc="varyhospitals")
# ## saves in files "out/varyhospitals-pra-FALSE-mN-trials100.Rdata"
table.efficiency.many.hospitals = function(many.m = c(4), n = 25, uniform.pra=F, trials=10) {
  for(m in many.m) 
    table.efficiency(m=m, 
                     sizes=c(n), 
                     uniform.pra=uniform.pra, 
                     trials=trials, 
                     filedesc="varyhospitals")
}

## Searches for files matching "varyhospitals-" and then 
## creates ONE png file 2x2 boxplot, where x-axis = #hospitals.
table.efficiency.many.to.graph = function() {
  files = list.files(path="out/", pattern="varyhospitals.*?-m", full.names=T)
  # boxplot to print (should be a list over hospital sizes)
  # i.e.  = [[mech.id]][[m.size]]
  boxplot.obj = list()
  
  for(filename in files) {
    x = filename
    m = regexec(text=x, pattern="-m(.*?)-")
    m.size = as.integer (  regmatches(x, m)[[1]][2]  )
    m.size.str = sprintf("%d", m.size)
    load(filename)
    # assume:  results. = [[mech]][[size]] = c(utility)
    for(mech.id in names(results) ) {
      
      if(length(results[[mech.id]]) != 1)
        stop("results should have size 1 (just one hospital size)")
      
      #loginfo(sprintf("Mechanism is %s m.size.str = %s", mech.id, m.size.str))
      if (!mech.id %in% names(boxplot.obj)) {
        boxplot.obj[[mech.id]] = list();
      }
      mech.utils = results[[mech.id]][[1]] # only one size.
      true.utils = results[["rCM_t"]][[1]]
      boxplot.obj[[mech.id]][[m.size.str]] = mech.utils / true.utils
    }
  }
  loginfo(names(boxplot.obj));
  png.file = sprintf("out/png/efficiency-varyhospitals.png")
  loginfo(sprintf("Saving in filename %s ", png.file))
  
  png(file= png.file)
  par(mfrow=c(2,2) )
  for(mech.id in names(boxplot.obj))
    boxplot(boxplot.obj[[mech.id]], ylim=c(0.6, 1.1), 
          main=sprintf("mech=%s", mech.id ), ylab="ratio (str/base)", xlab="#hospitals")
  dev.off();
}

















