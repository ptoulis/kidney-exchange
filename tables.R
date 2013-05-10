##  Tables of the paper. 
source("experiments.R")
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
    ##    alignment
    attachLine(sprintf("%s%s}","\\begin{tabular}{", 
                       paste(rep("c", length(columns)), collapse=" | ")))
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
##  Table 1.   Model assumptions.
##  n   |  #empirical matches |  theoretical matches
##   Synopsis.    D = table.assumptions(c(50,100,200, 300), trials=100)
##               table.to.tex(D, filename="out/assumptions.tex")
table.assumptions = function(sizes=c(50), trials=10) {
    ncols = 3 ## size-experimental matches - theoretical matches 
    pb = txtProgressBar(style=3)
    ## result matrix
    M = matrix(NA, nrow=length(sizes), ncol=ncols)
    ## standard error matrix
    SE =  matrix(NA, nrow=length(sizes), ncol=ncols)
    colnames(M) = c("n", "$\\mu_E(n)$", "$\\mu_T(n)$")
    
    ##  Total number of runs
    nruns = length(sizes) * trials
    
    for(i in 1:length(sizes)) {
        n = sizes[i]
        ##  x =  2 x trials
        x = sapply(1:trials, function(j) {
            val = (i-1) * trials + j
            setTxtProgressBar(pb, val/nruns)
            
            # Sample the RKE, then match
            rke = rrke(n)
            m = max.matching(rke)
            ## Empirical matches.
            empirical.matches = m$matching$utility
            theoretical.matches = 0.556 * n-0.338 * sqrt(n)-2
            return(c(empirical.matches, theoretical.matches))
        });
        
        M[i,] = c(n,  t(apply(x, 1, mean)) )
        SE[i,] = c(0, t(apply(x, 1, bootstrap) ))
    }   
    # text = 
    D = add.se(M, SE)
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

##   Table 3.   rCM and xCM scenarios. Do it size-by-size
table.mech = function(sizes, mech="rCM", sims=10) {
   ## Return matrix.
    k = 3
    M = matrix(NA, nrow=k * 4, ncol=1+length(sizes))
    colnames(M) = c("Hospital/size(n)", sizes)
    ## Scenarios
    dev = list(A=c(), B=c(1), C=c(2,3), D=c(1,2,3))
    
    ## Converts an intermed result to a vector.
    to.vector = function(x) {
        ##  x= [[hospital]]  = Data.frame
        v = c()
        for(i in 1:k) {
            b = x[[i]]
            b=  b[,5]
            index = which(b<=0)
            if(length(index)>0) {
                warning(sprintf("Found %d zeroes. LP output unstable.", length(index)))
            }
            ## Get only those who did not time out.
            b = b[b>0]
            
            ##  Get the "MatchTotal"
            ## Remove 0's
            
            m = mean(b)
            se = bootstrap(b)
            v[i] = sprintf("%.2f (%.2f)", m, se)
        }
        return(v)
    }
    
    ## For all scenarios
    for(j in 1:length(names(dev))) {
        scenario = names(dev)[j]
        deviating = c(dev[[scenario]])
        ## For all sizes
        for(i in 1:length(sizes) ) {
            size = sizes[i]
            
            x = e34.many.runs(k=k, sizes=c(size),mechanism=mech, 
                              deviating=deviating, 
                              sims=sims)
            # x[[size]][[hospital]]
            start.row=1+(j-1)*k
            end.row = j * k
            M[start.row:end.row, 1] = sapply(1:k, function(j) sprintf("H%d%s",j,
                                                                      ifelse(j %in% deviating,"d","")))
        
           M[start.row:end.row, i+1] = to.vector(x[[1]])
        }
    }
    return(M)
}



