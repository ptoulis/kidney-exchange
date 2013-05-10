# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges
source("rke2.R")
source("mechanisms.R")

## How many patients per ns step
NS.STEP = 50

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
## k = hospitals, n=size
## Returns  k x 5 matrix
##          ()
e34.one.run <- function(k, n, deviating, mechanism, 
                        verbose=F, 
                        use.rke.set=list()) {
    if(verbose)
        print(sprintf("E3: k=%d hospitals, %d deviating, %d truthful,  n=%d pairs", 
                      k,length(deviating), k-length(deviating),n))
    
    hospitals.deviating = rep(0,k)
    hospitals.deviating[deviating] = 1
    
    ##
    rke.set = NA
    if(length(use.rke.set)==0)
        rke.set = rrke.many(k,n)
    else rke.set = use.rke.set
    ## 
   
    ##  Maximum internal matches for every Hi
    max.internal.matches = rep(0, k)
    real.internal.matches = rep(0, k)
    
    ###     Deviating hospitals
    for(h in 1:k) {
        rke.h = rke.set[[h]]
        m = max.matching(rke.h)
        mids = m$matching$matched.ids
        remainder = remove.pairs(rke.h, mids)
        max.internal.matches[h] = length(mids)
        ## If i deviates then replace the report to the mechanism
        if(h %in% deviating) {
            rke.set[[h]] = remainder
            real.internal.matches[h] = max.internal.matches[h] 
        }
    }
    ## Pool all RKE's
    init = rep(NA,k)
    ##   k x 1   : From mechanism
    if(mechanism=="rCM")
        HospitalUtility = rCM(rke.set)
    else
        HospitalUtility = xCM(rke.set)
    #print(HospitalUtility)
    #print(hospitals.deviating)
    matched.vector = as.vector(HospitalUtility ) + real.internal.matches
    return(cbind(HospDev=hospitals.deviating,
                 MaxInternal=max.internal.matches, 
                 MatchInternal=real.internal.matches, 
                 MatchMech=as.vector(HospitalUtility ),
                 MatchTotal=matched.vector,
                GlobalMatch = rep(sum(matched.vector), k)))
    
}
## Many runs.
## Returns   x[[size]][[hospital]] =  sims x (N) matrix => N= (Dev, MaxInternal,...MatchTotal)
e34.many.runs <- function(k=3, sizes=c(50), 
                          sims=20, deviating=c(), mechanism="rCM") {
    col.names= colnames(e34.one.run(k, 20, deviating, mechanism))
    
    pb = txtProgressBar(style=3)
    results = list()
    progress.val = 0
    n.total = length(sizes) * sims
    
    for(i in 1: length(sizes)) {
        size = sizes[i]
        results[[i]] = list()
        for(j in 1:sims) {
            mat = e34.one.run(k,size, deviating=deviating, mechanism=mechanism, verbose=F)
            for(hid in 1:k) {
                if(length(results[[i]]) < hid)
                    results[[i]][[hid]] =  matrix(0, nrow=0,  ncol=length(col.names))
                colnames(results[[i]][[hid]]) = col.names
                results[[i]][[hid]] = rbind(results[[i]][[hid]], as.vector(mat[hid,]))
            }
            
            progress.val = progress.val+1
            setTxtProgressBar(pb=pb,value= progress.val/ n.total    )
        }
    }
    
    return(results)
}

## k = # hospitals
## ns.out = # of different hospital sizes
## sims  = # simulations
## mechanism = {rCM, xCM}
##    Implements experiments 3 and 4.
e34.scenario <- function(k, ns.out, sims, mechanism, scenario) {
    ## (a)
    if(scenario=="A")
        return(e34.many.runs(k, ns.out, sims, deviating=c(), mechanism= mechanism))
    
    ## (b)   One deviates
    if(scenario=="B")
        return( e34.many.runs(k, ns.out, sims, deviating=c(1), mechanism= mechanism) )
    
    ## (c)  All but one deviating
    if(scenario=="C")
        return(e34.many.runs(k,ns.out, sims, deviating=2:k, mechanism= mechanism))
    
    #  (D)  All deviate
    if(scenario=="D")
        return(e34.many.runs(k,ns.out, sims, deviating=1:k, mechanism= mechanism))
    
}



#######      EXPERIMENTS start HERE   ########


##  Basic table for  mu(n) = E M(n)  = expected # of matches
##
##   Return: out,    e.g. out[[1]] = list(avg, sd)
E1 <- function(ns.out=5, sims=100) {
    out= list()
    pb = txtProgressBar(style=3)
    val = 0
    total.N = ns.out * sims
    for(i in  1:ns.out) {
        n = NS.STEP * i
        
        reps = replicate(sims, {  rke = rrke(n)
                                  
                                  m = max.matching(rke)
                                  length(m$matching$matched.ids)});
        val = val+sims
        setTxtProgressBar(pb, value=val/total.N)
        out[[i]] = list(avg=mean(reps), sd=sd(reps))
    }
    close(pb)
    return(out)
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

###    Experiment 3 :   Incentives and efficiency for rCM 
##   results[[i]][[j]]  =  sims x K  matrix

##   i = hospital size (20 * i)    j=hospital id
##   sims = no. of simulations
##   K = data columns = 
##       (HospitalDeviating,   MaxInternal,  MatchInternal,  MatchMechanism, MatchTotal)
##  1.  To get the mean matches of hospital 3 of size 40 write:
##            apply(results[[2]][[3]], 2, mean)
##  2.  To compute the bootstrap distribution:
##            v = res[[2]][[3]]
##            reps = replicate(1000, { mean(sample(v, 1000, replace=T))})
##            hist(reps)
E3 <- function(k=3,ns.out=3,sims=50, scenario="A") {
    
    results = e34.scenario(k, ns.out, sims, "rCM", scenario)
    return(results)
}


## Experiments on xCM  (See E3() )
E4 <- function(k=3,ns.out=3,sims=50, scenario="A") {
    
    results = e34.scenario(k, ns.out, sims, "xCM", scenario)
    return(results)
}

#
###
###           LaTeX   functions    ########
#
#e34.one.run <- function(k, n, deviating, mechanism, 
#                        use.rke.set=list()) {
TripleExperiment <- function(k=3, n=300, trials) {
    results = list()
    pb = txtProgressBar(style=3)
    
    for(i in 1:trials) {
        rke.set = rrke.many(k, n)
        
        results[[i]] = list()
        out1 = e34.one.run(k,n, c(), "rCM", use.rke.set = rke.set)
        out2 = e34.one.run(k, n, c(), "xCM", use.rke.set = rke.set)
        out3 = e34.one.run(k,n, 1:3, "xCM", use.rke.set = rke.set)
        
        if(out2[1,6] > out1[1,6])
        {
            print("oops")
            print(out1)
            print(out2)
            save(rke.set, file="out/problematic-triple.Rdata")
            readline()
            
        }
        
        results[[i]][["rCMa"]]= out1
        results[[i]][["xCMa"]]= out2
        results[[i]][["xCMd"]]= out3
        
        
        save(results, file="out/triple-experiment.Rdata")
        setTxtProgressBar(pb, value=i/trials)
    }
    return(results)
}

### From E1 obejct to LaTeX
e1.to.latex<- function(y) {
    ns = length(y) ## sizes
    lines = c()
    fileConn = file("out/output.txt",open="a")
    
    attachLine <- function(s) {
        lines[length(lines)+1] <<- s
    }
    attachLine("")
    attachLine("------")
    attachLine("\\begin{center}")
    columns = c("n", "\\#matches", "$\\mu(n)$")
    ##    alignment
    attachLine(sprintf("%s%s}","\\begin{tabular}{", paste(rep("l", length(columns)), collapse=" | ")))
    tHeader = columns
    attachLine(sprintf("%s\\\\ \\hline", paste(tHeader, collapse=" & ") ))
    
    
    for(i in  (1:ns)) {
        n = NS.STEP * i
        obj = y[[i]]
        val = sprintf("%.2f (%.2f)", obj$avg, obj$sd)
        mu.n = sprintf("%.2f", 0.556 * n - 0.338 * sqrt(n) - 2)
        values = c(n, val, mu.n)
        attachLine(sprintf("%s \\\\",paste(values, collapse=" & ")))
        
    }
    attachLine("\\end{tabular}")
    attachLine("\\end{center}")
    writeLines(lines,con=fileConn)
    close(fileConn)
    
    }
## Convert experiments output to LaTeX
e2.to.latex <- function(y , 
                        suffix="default",
                        filename) {
  
    ns = length(y) ## sizes
    sample.run = e2.once(size=20)
    
    lines = c()
    fileConn = file(filename,open="w")
    
    attachLine <- function(s) {
        lines[length(lines)+1] <<- s
    }
    attachLine("")
    attachLine("---  E2  ---")
  
    attachLine("\\begin{table}[ht]")
    attachLine("\\centering")
    columns = colnames(sample.run)
    ##    alignment
    attachLine(sprintf("%s%s}","\\begin{tabular}{", paste(rep("c", length(columns)), collapse=" | ")))
    tHeader = c(columns)
    attachLine(sprintf("%s\\\\ \\hline", paste(tHeader, collapse=" & ") ))
    
    
    for(i in  (1:ns)) {
        experiment = y[[i]]
        n = experiment[1,1] ## n= size of hospital
        sims = nrow(experiment)
        values = c()
        ##  y[[i]] = data frame of violations
        std.errors = as.vector(apply(y[[i]],2, function(x) sd(x)/sqrt(sims)))
        means = as.vector(apply(y[[i]], 2, mean))
        new.values = sapply(1:length(means), function(x)
                            sprintf("%.2f (%.2f)", means[x], std.errors[x]))
        
        values = c(values, new.values)
        attachLine(sprintf("%s \\\\",paste(values, collapse=" & ")))
        
    }
    attachLine("\\end{tabular}")
    attachLine("\\vspace{0.4cm}")
    attachLine("\\caption{Caption.}")
    attachLine("\\end{table}")
    writeLines(lines,con=fileConn)
    close(fileConn)
    }
## Convert experiments output to LaTeX
##  y  = output from #3 or #4 experiment.
e34.to.latex <- function(y, suffix="default", tableDesc="E34") {
    
    ns = length(y) ## sizes
    n.hospitals = length(y[[1]])
    lines = c()
    fileConn = file("out/output.txt", open="a")
    getHospitalName <- function(hid) {
        m = as.data.frame(y[[1]][[hid]])
        if(m$HospDev[1]==1)
            return(sprintf("H-%d$_d$", hid))
        return(sprintf("H-%d", hid))
    }
    attachLine <- function(s) {
        lines[length(lines)+1] <<- s
    }
    attachLine("")
    attachLine(sprintf("--  %s ---", tableDesc))
    attachLine("\\begin{center}")
    attachLine(sprintf("%s%s}","\\begin{tabular}{", paste(rep("l |", n.hospitals+1), collapse=" ")))
    tHeader = c("n", sapply(1:n.hospitals, function(i) getHospitalName(i)))
    attachLine(sprintf("%s\\\\ \\hline", paste(tHeader, collapse=" & ") ))
    
    
    for(i in  (1:ns)) {
        n = NS.STEP * i
        values = c(n)
        for(h in 1:n.hospitals) {
            m = as.data.frame(y[[i]][[h]])
            values[length(values)+1] =
                sprintf("%.2f (%.2f)", mean(m$MatchTotal), sd(m$MatchTotal))
        }
        attachLine(sprintf("%s \\\\",paste(values, collapse=" & ")))
        
    }
    attachLine("\\end{tabular}")
    attachLine("\\end{center}")
    writeLines(lines,con=fileConn)
    close(fileConn)
}
## old.



