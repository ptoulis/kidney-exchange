# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges

##  pair = { donor, patient, pra, type }     

Blood.Types  <- c("O", "A", "B", "AB")
Blood.Codes  <- c(1, 2, 3, 6)
PRA.vals    <- c(0.2)
## Given a pair -> [ compatible ABO]
possible.matches <-function(pair) {
    donor = pair$donor
    patient = pair$patient
    dons = sapply(Blood.Codes, function(i) as.blood.code(patient)%%i==0)
    dons = Blood.Types[which(dons==1)]
    pats = sapply(Blood.Codes, function(i) i%%as.blood.code(donor)==0)
    pats = Blood.Types[which(pats==1)]
    all.types=c()
    for(i in dons)
        for(j in pats)
            all.types <-  c( all.types, sprintf("%s-%s", i,j))
    return(all.types)     
}
#   From type (e.g. "O") to numeric code
as.blood.code <- function(b.type) {
    if(is.numeric(b.type)) return(b.type)
    if(! is.vector(b.type))
       return(Blood.Codes[ which(Blood.Types==b.type)])
    
    return( sapply(1:length(b.type), function(j) Blood.Codes[ which(Blood.Types==b.type[j])]) )
}
as.blood.type <- function(b.code) {
    i = which(Blood.Codes==b.code)
    return(Blood.Types[i])
}

#  Sample blood type.
rblood <- function(n, probs=c(50,30,15,5)) {
    return( sample( Blood.Types, size=n, replace=T, prob=probs))
}

## Sample a PRA value. For now this is constant to 0.2
rpra <- function(n, is.uniform=F) {
    if(is.uniform)
        return( sample(PRA.vals, n, replace=T))
    warning("Using non-uniform PRA values.")
    non.uniform.vals = c(0.05, 0.45, 0.9)
    return(sample(non.uniform.vals, size=n, replace=T, prob=c(0.7, 0.2, 0.1)))
}
## Sample pairs.
rpairs <- function(n) {
    p = list(donor=rblood(n), 
             patient=rblood(n),
             pra = rpra(n))
    p$type = pair.type(p)
    return(p)
}

## The type as overdemanded (O), under-demanded (U), reciprocal (R) and 
## self-demanded (S)
pair.type <- function(pair) {
    bp = as.blood.code(pair$patient)
    bd = as.blood.code(pair$donor)
    b = (bp %% bd == 0)
    if(b) {
        if(bp==bd) return("S")
        return("O")
    } else {
        if(bp==2 && bd==3) return("R")
        if(bp==3 && bd==2) return("R")
        return("U")
    }
}

## Returns the color of the pair
## Used for plotting.
pair.color <- function(pair) {
    type= pair.type(pair)
    if(type=="U") return("gray")
    if(type=="O") return("green")
    if(type=="R") return("yellow")
    if(type=="S") return("cyan")
}

pair.self.compatible <- function(pair) {
    ## 1. Donor can donate to patient?
    b = (as.blood.code(pair$patient) %% as.blood.code(pair$donor) == 0)
    ## 2. Crossmatch avoided?
    not.cross = runif(1) >  pair$pra
    ## return TRUE if blood type compatible AND no crossmatch occured
    return(b && not.cross)
}
pairs.compatible <- function(pair1, pair2) {
    not.cross = (runif(1) > pair1$pra) && (runif(1) > pair2$pra)
    if(! not.cross) return(F)
    ## Check if cross blood-type compatible
    b1 = (as.blood.code(pair1$patient) %% as.blood.code(pair2$donor) == 0)
    if(! b1) return(F)
    b2 = (as.blood.code(pair2$patient) %% as.blood.code(pair1$donor) == 0)
    if(! b2) return(F)
    ## Check if no crossmatch occurs
    return(T)
}

## Pair code is from 1 to 16
decode = function(pair.code) {
    i = pair.code
    pat = ifelse(i %% 4==0, 4, i%%4); 
    donor = 1+(i-pat)/4;
    return(list(donor=donor, patient=pat))
}
compatible.patients <- function(donor.code) {
    x = c()
    if(donor.code==1) x = c(1,1,1,1)
    if(donor.code==2 )x = c(0,1,0,1)
    if(donor.code==3) x = c(0,0,1,1)
    if(donor.code==4) x = c(0,0,0,1)
    
    return(  rep(x, 4) )
}
compatible.donors <- function(pat.code) {
    x = c()
    if(pat.code==1) x = c(rep(1,4), rep(0, 12))
    if(pat.code==2) x = c(rep(1,8), rep(0, 8))
    if(pat.code==3) x = c(rep(1,4), rep(0,4), rep(1,4), rep(0,4))
    if(pat.code==4) x = rep(1,16)
    return(x)
}
##  returns the codes that are b-compatible with the pair.
compatible.codes <- function(pair1) {
    dec = decode(pair1)
    mask1 = compatible.patients(dec$donor)
    #print(mask1)
    
    mask2 = compatible.donors(dec$patient)
    #print(mask2)
    return( which(mask1 * mask2==1) )
}

##
##
##
##
###      Added: Jan. 30  --- New functions for sampling and pooling. Faster!
## 
##     pair.code = integer in {1,2.,,....16} for each pair O-O, O-A, ...
##     New object = {   compact => compact matrix,  pc => pair codes   }
##
##  Samples a binary matrix, based on probabilities of prob.matrix
rbin.matrix<- function(prob.matrix) {
    n = nrow(prob.matrix)
    a.vec = as.vector(prob.matrix)
    y = rbinom(length(a.vec), size=1, prob=a.vec)
    
    rm(a.vec)    
    ret.matrix = matrix(y, nrow=n)
    rm(y)
    #gc()
    
    return(ret.matrix)
}

##  Samples P = binary cross-match matrix
## Can have 1 or 2 hospitals. If same we need to treat the diagonal elements differently.
## If equal to 1, then they can self match.
sample.bin.pra.matrix = function(pra.probs1, pra.probs2, same.hospital=T,  verbose=F) {

    
    Us1 = 1-pra.probs1
    Us2 = 1-pra.probs2
    pras.matrix = Us1 %*% t(Us2)
    
    
    if(same.hospital)
        pras.matrix = pras.matrix - diag( Us1^2 )
    ## Why -diag?   Because if we are in the same hospital we know 
    ## that the pairs cannot be self-matched.

    ## M = pras matrix   Mij = (1-pi) * (1-pj) and Mii = (1-pi)
   
    ret.matrix = rbin.matrix(pras.matrix)
    rm(list=c("pras.matrix", "Us1", "Us2"))
    d = diag(ret.matrix)
    
    ##  Get the upper triangular and symmetrize.
    up.ret = upper.tri(ret.matrix) * ret.matrix
    ret.matrix = up.ret + t(up.ret) + diag(d)
    ###   Symmetrize
    if(verbose) {
        print("Length of binom calls ")
        print(length(pras.vec))
        print("PRA probs matrix")
        print(pras.matrix)
        print("PRA vector")
        print(pras.vec)
        print("Y binomial ")
        print(y)
        print("PRA compatibility matrix")
        print(ret.matrix)
    }
  
    rm(list=c("up.ret", "d"))

    return( ret.matrix)
}

###   Returns B = binary blood-type compatibility matrix
get.bin.blood.matrix <- function(pair.codes, verbose=F) {
    n = length(pair.codes)
    
    if(verbose) {
        print("Sampled pair codes")
        print(pair.codes)
    }
    ### Sample the binary blood-type compatibility matrix.
    pair.ids = 1:n
    B = matrix(0, nrow=n, ncol=n)
    cursor = 1
    dat = as.data.frame(table(pair.codes))
    if(verbose)
        print(dat)
    for(row.i in 1:nrow(dat)) {
        #  1.  get one pair code  Y, e.g. Y = 7  (from 1 to 16)
        p.code = as.numeric(as.vector(dat$pair.codes))[row.i]
        # 2.  index which are = Y
        i.range = which(pair.codes==p.code)
        # 3. Codes of nodes who are compatible with Y
        compat.vec = compatible.codes(p.code)
        # 4. index of sampled pair codes = compatibile ones
        j = which(pair.codes %in% compat.vec)
        
        if(verbose) {
            print(sprintf("Pair code = %d", p.code))
            print(i.range)
            print("Compatible vector")
            print(compat.vec)
            print("Compatible pairs")
            print(j)
            
        }
        # 5.  Set those edges equal to 1
        B[i.range, j] = 1
    }
    # cleanup
    rm(pair.codes)
    rm(pair.ids)
    return(B)
}

self.matched.pair = function(pair.code) {
    return(pair.code %in% compatible.codes(pair.code))
}

## Samples pairs and their PRAs
rpairs.new <- function(n, sample.pra.fn=rpra, verbose=F) {
    out.codes= c()
    out.pras= c()
    basic.p = c(0.5, 0.3, 0.15, 0.05)
    all.pairs = 1:16
    probs = sapply(all.pairs, function(i) { dec = decode(i); 
                                       basic.p[dec$patient] * basic.p[dec$donor]})
    
    while(length(out.codes)<n) {
        no.samples =  n+100
        pair.codes = sample(all.pairs, size= no.samples, replace=T, prob=probs)
        pras = sample.pra.fn( no.samples)
        no.cross = rbinom( no.samples, size=1, prob=1-pras)
        b = as.numeric( sapply(pair.codes, function(i) self.matched.pair(i)) )
        j = which( b * no.cross==0)
        if(length(j)>0)
         {
            out.codes = c(out.codes, pair.codes[j])
            out.pras = c(out.pras, pras[j])
        }
    }
    shuffled = sample(1:length(out.codes), size=n, replace=F)
    return( list(codes=out.codes[shuffled], pras=out.pras[shuffled]))
}


## Converts a pair code to a pair
pair.code.to.pair <- function(pair.code) {
    dec=  decode(pair.code)
    return(list(donor= Blood.Types[dec$donor], patient=Blood.Types[dec$patient]))
}
to.rke <- function(new.rke) {
    rke = list(pairs= list(), graph=graph() )
    # 1. Set the graph
    rke$graph = graph.adjacency(A=new.rke$A)
    
    # 2. Set the pairs.
    for(i in 1:length(new.rke$pair.codes)) {
        pair.code = new.rke$pair.codes[i]
        pair.pra = new.rke$pras[i]
        
        i = length(rke$pairs)+1
        pair = pair.code.to.pair(pair.code)
        rke$pairs[[i]] = list(donor=pair$donor, 
                              patient=pair$patient, 
                              pra=pair.pra,
                              type=pair.type(pair))
        
    }
    return(rke)
}




