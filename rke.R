# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges
## Jan 2013, new version of rke.R
rm(list=ls())
source("terminology.R")

rrke.many <- function(m=3, n=60, uniform.pra, verbose=F) {
  # Create an RKE objects with many hospitals.
  #
  # Args: m = #hospitals, n = #pairs
  # Returns: RKE object
  if (m == 0 | n == 0)
    return(empty.rke())
  all.pairs = 1:(m*n)
  get.pair.ids <- function(hid) (1+(hid-1)*n):(hid*n)
  rke.all = rrke(n, pair.ids=get.pair.ids(1), uniform.pra=uniform.pra, hospital.id=1)
  if (m >= 2)
    for(hospital.id in 2:m) {
      rke = rrke(n, pair.ids=get.pair.ids(hospital.id),
                 uniform.pra=uniform.pra, hospital.id=hospital.id)
      rke.all = append.rke(rke.all, rke, verbose=verbose)
    }
  return(rke.all)
}

## Total number of pairs in the graph
rke.size <- function(rke) nrow(rke$pairs)
rke.remove.pairs <- function(rke, pair.ids) {
  ##  Given an RKE  (1) subtract (2) return remainder.
  # Used to represent deviation strategies ("hide")
  CHECK_MEMBER(pair.ids, rke.pair.ids(rke))
  rke$pairs = subset(rke$pairs, !is.element(pair.id, pair.ids))
  rke$edges = sample.pairs.edges(pairs=rke$pairs, respect.edges=rke$edges)
  rke$A = map.edges.A(rke$edges)
  return (rke)
}

rke.add <- function(rke1, rke2, verbose=F) {
  # Combines two RKE objects into one.
  warning("Append is not unit-tested")
  CHECK_rke(rke1)
  CHECK_rke(rke2)
  if (nrow(rke1$pairs) == 0)
    return(rke2)
  if (nrow(rke2$pairs) == 0)
    return(rke1)
  CHECK_TRUE(all(!is.element(rke1$pairs$pair.id, rke2$pairs$pair.id)),
             msg="RKEs should not have the same pairs.")
  rke.all = empty.rke()
  all.edges = rbind(rke1$edges, rke2$edges)
  rke.all$pairs = rbind(rke1$pairs, rke2$pairs)
  rke.all <- update.rke.new.pairs(rke=rke.all, keep.edges=all.edges)
  return(rke.all)
}

rke.keep.pairs = function(rke, pair.ids) {
  # Keeps the specified ids and removes all others.
  # Returns a <rke> object.
  all.pairs = rke.pair.ids(rke)
  rm.pairs = setdiff(all.pairs, pair.ids)
  return(rke.remove.pairs(rke, rm.pairs))
}

rke.edge.ids = function(rke) rke$edges$edge.id
rke.pair.ids = function(rke) rke$pairs$pair.id

rke.2way.cycles <- function(rke) {
  # Returns a K x 2 matrix of 2-way cycles in the rke object.
  # (K = #2-way cycles)
  A = rke$A
  pair.ids = as.numeric(rownames(A))
  cycles.two = apply(which((A * t(A)) == 1, arr.ind=T), 2, function(i) pair.ids[i])
  x = apply(cycles.two, 1, function(r) r[2] > r[1])
  return(cycles.two[which(x), ])
}

rke.3way.cycles <- function(rke) {
  # Given an adjacency matrix, it will return the 3-way cycles.
  #
  # Args: A adjacency matrix
  #
  # Returns: k x 3 matrix of cycles. Each row is (i,j,k) indicates i->j->k->i cycle
  A = rke$A
  cycles.3way = matrix(0, nrow=0, ncol=3)
  num.nodes = nrow(A)
  subset.nodes = 1:num.nodes
  have.inlinks = which(colSums(A) >= 1)
  have.outlinks = which(rowSums(A) >= 1)
  subset.nodes = intersect(have.inlinks, have.outlinks)
  for (i in subset.nodes) {
    Ni = which(A[, i] == 1)  # ids of inlinks
    if (length(Ni) == 0)
      stop(sprintf("Node %d should have inlinks.", i))
    BiNi = matrix(0 , nrow=num.nodes, ncol=length(Ni))
    for (j in 1:length(Ni))
      BiNi[Ni[j], j] <- 1
    rm(j)
    # BiNi = binary representation of inlinks of i: n x Ni (Ni=#inlinks)
    # [, j] = (0,0, ...1, 0,0), 1 at j-th inlink
    Oi = which(A[i, ] == 1)  # outlinks of i
    if (length(Oi) == 0)
      stop(sprintf("Node %d should have outlinks", i))
    Oi.outlinks = A[Oi, ]  # outlinks of outlinks  (Oi x n) matrix.
    OiBi = Oi.outlinks %*% BiNi
    # OB = oulinks-binary-inlinks matrix: 
    # OB_jk = 1 iff there is a 3-way cycle i-> Oj ->Ik
    # i.e. between i, its j-th outlink and its k-th inlink.
    found.cycles = which(OiBi == 1, arr.ind=T)
    if (length(found.cycles) > 0)
      for (j in 1:nrow(found.cycles)) {
        row = found.cycles[j, 1]
        col = found.cycles[j, 2]
        cycles.3way <- rbind(cycles.3way, c(i, Oi[row], Ni[col]))
      }
  }
  # Remove the permutations
  # Take only those cycles where the ids only increase or decrease
  # i.e. 3,2,1  or 1,2,3
  # In 3-way cycles one will always exist
  x = t(apply(cycles.3way, 1, diff))
  x = apply(x, 1, function(r) all(r < 0) || all(r > 0))
  ret = cycles.3way[which(x), ]
  # map back to pair ids.
  pair.ids = as.numeric(rownames(A))
  return (apply(ret, 2, function(i) pair.ids[i]))
}

plot.rke = function(rke, vertex.size=20) {
  library(igraph)
  g = graph.adjacency(rke$A, mode="directed")
  V(g)$color = rke$pairs$pair.color
  V(g)$label = str_c("H", rke$pairs$hospital, "#", rke$pairs$desc, rke$pairs$pair.id)
  par(mar=c(0,0,0,0))
  plot.igraph(g,layout=layout.auto, vertex.size=vertex.size)
}