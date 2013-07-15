# matching with 3-way exchanges.

A = matrix(0, nrow=6, ncol=6)
A[1, 3] <- 1
A[2, c(1,5)] <- 1
A[3, c(2,4,5)] <- 1
A[4, 3] <- 1
A[5, 3] <- 1

library(igraph)
g = graph.adjacency(A, mode="directed")

get.3way.cycles <- function(A) {
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
  return(cycles.3way[which(x), ])
}

get.2way.cycles <- function(A) {
  cycles.two = which((A * t(A)) == 1, arr.ind=T)
  x = apply(cycles.two, 1, function(r) r[2] > r[1])
  return(cycles.two[which(x), ])
}