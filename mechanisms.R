## Mechanisms


## Mechanisms return   k x 1  matrix  with the matches for every hospital.
get.hospitals.utility <- function(rke.all, m.all) {
    k = length(unique(rke.all$hospital))
    vals = matrix(0, nrow=k, ncol=1)
    matched.ids = m.all$matching$matched.ids
    matched.hospitals =rke.all$hospital[matched.ids]
    for(hid in 1:k) {
       vals[hid,1] = length(which(matched.hospitals==hid)) 
    }
    return(vals)
}

## Implementation of rCM

rCM <- function(rke.list) {
    k = length(rke.list)
    HospitalUtility = matrix(0, nrow=k, ncol=1)
    rke.all = pool.rke(rke.list=rke.list)
    ##  Set a time limit?
    m.all =  max.matching(rke.all)
   
    HospitalUtility = get.hospitals.utility(rke.all, m.all)
    return(HospitalUtility)
}

## Implementation of practical xCM
## 
xCM <- function(rke.list) {
    original.list = rke.list
    k = length(rke.list)
    IR.cons= matrix(0, nrow=k, ncol=2)
    
    ## Utility = K x 1 {  total matched }
    HospitalUtility = matrix(0, nrow=k, ncol=1)
    to.use.cutoff = F # k * get.size(rke.list[[1]]) > 1000
    
    for(hid in 1:k) {
        rke.h = rke.list[[hid]]
        A = get.model.A(rke.h)
      
        m.h = max.matching(rke.h, 
                           regular.matching=T, 
                           use.cutoff = to.use.cutoff)
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
    ###  Warning. For some, this takes very long time. We need to bypass.
#     m.all = evalWithTimeout(max.matching(rke.all,IR.constraints=IR.cons), 
#                             timeout=10, 
#                             onTimeout="silent")
#     if(length(m.all)==0)
#        return(NULL)
    m.all = max.matching(rke.all,  IR.constraints=IR.cons)
    
    HospitalUtility = HospitalUtility + get.hospitals.utility(rke.all, m.all)
                                 
    
    return(HospitalUtility)
    
}