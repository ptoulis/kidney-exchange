##  Tables of the paper. 
## Code to save into latex table.
library(xtable)

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

###      FUNCTIONS TO CREATE THE TABLES
##  Table 1.   Theorem 1 (expected value of maximum matchings.)
##    n | selifh | together | m(n)   | surplus
##   Synopsis.    D = table.matchings(sizes = c(50,100,200, 300), m=2, trials=100)
##               table.to.tex(D, filename="out/assumptions.tex")
table.matchings = function(sizes=c(50), 
                           m= 5, 
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
        n = sizes[i]
        ##  x =  2 x trials
        x = sapply(1:trials, function(j) {
            val = (i-1) * trials + j
            setTxtProgressBar(pb, val/nruns)
            
            # 0. Sample rke's 
            rke.list = rrke.many(m=m, n=n, uniform.pra=T)
            
            # 1.   Selfish (each one matches in isolation) 
            selfish.matches = c()
            for(hid in 1:length(rke.list)) {
              selfish.matches[hid] = max.matching(rke.list[[hid]])$matching$utility; 
            }
            selfish.matches = sum(selfish.matches)
            
            # 2. Pool and match together
            rke.all = pool.rke(rke.list)
            m.all = max.matching(rke.all)
            together.matches = m.all$matching$utility
            
            # 3. Theoretical matches
            together.theoretical = 0.556 * n *m -0.338 * sqrt(n * m)- 2
            
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
##  Call   table.violations(c(50, 100..),  sims=100)
table.violations = function(sizes=c(50), sims=10) {
     x = E2(sizes=sizes, sims=sims)
     ##   x[[size]]
     ##  Results matrix.
     M = matrix(NA, nrow=length(sizes), ncol=7)
     SE = matrix(NA, nrow=length(sizes), ncol=7)
     for(i in 1:length(sizes)) {
         result = x[[i]]
         
         M[i,] = apply(result, 2, mean)
         SE[i,] = apply(result, 2, bootstrap)
     }
     D = combine.mat.se(M, SE)
     return(D)
}

## Compute the relative gain for _hid.interest_ 
# mechanism _mech_ under _setup_, for two different strategies
# Returns:  2 x trials   matrix of utiliies.
relative.gain = function(kpd1, kpd2, mech, hid)
{
  U1 = Run.Mechanism(kpd1, mech=mech)
  U2 = Run.Mechanism(kpd2, mech=mech)
  c(U1[hid], U2[hid])
}

##  Run a scenarion for a specific mechanism.
relative.gain.scenario = function(scenario, mech, m, n, trials) {
  
  h1.str.list = list("rCM" = list(A="c", B="c", C="c", D="c"),
                     "xCM" = list(A="c", B="r", C="c", D="r"),
                     "Bonus"=list(A="c", B="r", C="c", D="r") )
  
  others.list = list("rCM" = list(A="t", B="c", C="t", D="c"),
                "xCM" = list(A="t", B="t", C="t", D="t"),
                "Bonus"=list(A="t", B="t", C="t", D="t"))
  
  uniform.pra.list = list(A=T, B=T, C=F, D=F)
  
  ## define the parameters based on mechanism + scenario
  h1.strategy = h1.str.list[[mech]][[scenario]]
  others = rep(others.list[[mech]][[scenario]],  m-1)
  uniform.pra = uniform.pra.list[[scenario]]  
  
  A = matrix(NA, nrow=2, ncol=trials)
  
  for(i in 1:trials) {
    rke.list = rrke.many(m=m, n=n, uniform.pra=uniform.pra)
    rke.all = pool.rke(rke.list)
    ## e.g. Scenario A ->  str1 = "ctt"        str2 = "ttt"
    str1 = paste(c(h1.strategy, others), collapse="")
    ## str2 =  baseline strategy profile.
    str2 = paste(c("t", others), collapse="")
    
    kpd1 = kpd.create(rke.list, rke.all=rke.all, str1)
    kpd2 = kpd.create(rke.list, rke.all=rke.all, str2)
    
    utils = relative.gain(kpd1, kpd2, mech=mech, hid=1)
    A[1,i] = utils[1]
    A[2,i] = utils[2]
  }
  
  return(A)
}

##  Create tables 3,4,5 : Mechanisms.
table.mechs = function(mech, m=3, sizes=c(20), trials=10) {
  ## Return matrix.
  results = list()  
  scenarios = c("A", "B", "C", "D")
  N = length(scenarios) * length(sizes)
  pb = txtProgressBar(style=3)
  cnt = 0
  for(scen in scenarios)   {
    results[[scen]] = list()
    for(n in sizes) {
      
      results[[scen]][[sprintf("%d",n)]] = relative.gain.scenario(scenario=scen,
                                                                  mech=mech,
                                                                  m=m, n=n, trials=trials)
      cnt = cnt + 1
      setTxtProgressBar(pb, value=cnt / N)
    }
  }
  
  return(results)
}



