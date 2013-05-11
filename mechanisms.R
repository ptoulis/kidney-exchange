##   Panos Toulis  ptoulis@fas.harvard.edu
## Contains code implementing different KPD mechanisms. 
##  Currently supported:   rCM,  xCM
## A mechanism receives a combined donor-patient graph and outputs a matching. 
## (note: code returns a vector of utilities (matches/hospital)


## Given the matching, and the input graphs, 
##  calculate the utilities/hospital 
##  Returns    k x 1  vector,     with the matches for every hospital.
get.hospitals.utility <- function(rke.all, m.all) {
    ## find no. of hospitals 
    m = length(unique(rke.all$hospital))
    vals = matrix(0, nrow=m, ncol=1)
    ## get the pair ids in the match.
    matched.ids = m.all$matching$matched.ids
    ## get the hospitals for every pair.
    matched.hospitals =rke.all$hospital[matched.ids]
    
    ## Count matches/hospital. TO-DO(ptoulis): Use "table" here.
    for(hid in 1:m) {
       vals[hid,1] = length(which(matched.hospitals==hid)) 
    }
    return(vals)
}

##  Given an RKE and a type of deviation, report back the pairs which are matched internally.
play.strategy <- function(rke, type="truthful") {
  if(type=="truthful")
    return(c())
  if(type=="canonical") {
    ## TO-DO(ptoulis): If the max matching is timed-out, this will return empty match
    ## which is equivalent to being truthful. Needs a fix?
    m = max.matching(rke)
    return(m$matching$matched.ids)
  }
  stop("Not Implemented")
}

## Implementation of rCM
## deviating: list of hospitals which deviate.
##
##  Return:  kx1   matrix of utilities
rCM <- function(rke.list, strategies) {
    m = length(rke.list)
    HospitalUtility = matrix(0, nrow=m, ncol=1)
    
    
    for(hid in 1:m) {
      rke.h = rke.list[[hid]]
      matched.internally = play.strategy(rke.h, type=strategies[hid])
      HospitalUtility[hid,1] = length(matched.internally)
      rke.list[[hid]] = remove.pairs(rke.h, matched.internally)
    }
    
    ## 0.  Pool all the reports.
    rke.all = pool.rke(rke.list=rke.list)
    
    ## 1. Simply calculate a maximum-matching (this will shuffle the edges by default)
    m.all =  max.matching(rke.all)
    
    # 2. Compute the utility.
    HospitalUtility = HospitalUtility + get.hospitals.utility(rke.all, m.all)
    
    return(HospitalUtility)
}

## Implementation of practical xCM
## 
xCM <- function(rke.list, deviating=c()) {
    original.list = rke.list
    m = length(rke.list)
    IR.cons= matrix(0, nrow=m, ncol=2)
    
    ## Utility = K x 1 {  total matched }
    HospitalUtility = matrix(0, nrow=m, ncol=1)
    
    for(hid in 1:m) {
        rke.h = rke.list[[hid]]
        A = get.model.A(rke.h)
      
        m.h = max.matching(rke.h, 
                           regular.matching=T)
        # The edges containing O's will be removed in the matching
        mEdges = m.h$matching$matched.edges
        mOedges = intersect(filter.edges.by.type(rke.h,  "O","*"), 
                               mEdges) 
        mRRedges = intersect(filter.edges.by.type(rke.h,"R","R"),
                                mEdges)
        mSSedges = intersect(filter.edges.by.type(rke.h, "S","S"),
                                mEdges)
        
        
        
        remove.pairs = get.incident.nodes(A, mOedges)
        
        ## Remove all matched O pairs (and their matched) from the RKE
        rke.list[[hid]] = remove.pairs(rke.h, remove.pairs)
        
        IR.cons[hid,1] = 2 * length(mRRedges)
        IR.cons[hid, 2] = 2 * length(mSSedges)
        
        HospitalUtility[hid,1] = length(remove.pairs)
    }
    
    rke.all = pool.rke(rke.list=rke.list)
    m.all = max.matching(rke.all,  IR.constraints=IR.cons)
    
    HospitalUtility = HospitalUtility + get.hospitals.utility(rke.all, m.all)

    return(HospitalUtility)
}