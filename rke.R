# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges
source("terminology.R")
library(plyr)

rrke.pool <- function(m=3, n=60, uniform.pra, verbose=F) {
  # Create a "rke pool" object.
  #
  # Args: m = #hospitals, n = #pairs
  if (m == 0 | n == 0)
    return(empty.rke())
  all.pairs = 1:(m*n)
  get.pair.ids <- function(hid) (1+(hid-1)*n):(hid*n)
  rke.all = rrke(n, pair.ids=get.pair.ids(1), uniform.pra=uniform.pra, hospital.id=1)
  rke.list = list()
  rke.list[[1]] = rke.all
  if (m >= 2)
    for(hospital.id in 2:m) {
      rke = rrke(n, pair.ids=get.pair.ids(hospital.id),
                 uniform.pra=uniform.pra, hospital.id=hospital.id)
      rke.all = rke.add(rke.all, rke, verbose=verbose)
      rke.list[[hospital.id]] <- rke
    }
  return(list(rke.list=rke.list, rke.all=rke.all))
}

rke.equal <- function(rke1, rke2) {
  all(table(rke1$pairs$desc) == table(rke2$pairs$desc))
}

## Total number of pairs in the graph
rke.size <- function(rke) {
  CHECK_EQ(nrow(rke$pairs), length(unique(rke$pairs$pair.id)))
  nrow(rke$pairs)
}

rke.remove.pairs <- function(rke, rm.pair.ids) {
  ##  Given an RKE  (1) subtract (2) return remainder.
  # Used to represent deviation strategies ("hide")
  CHECK_MEMBER(rm.pair.ids, rke.pair.ids(rke))
  if (length(rm.pair.ids) == 0)
    return (rke)
  rke$pairs = subset(rke$pairs, !is.element(pair.id, rm.pair.ids))
  rke = rke.update.new.pairs(rke=rke, keep.edges=rke$edges)
  rke$A = map.edges.adjacency(rke$edges, rke$pairs$pair.id)
  return (rke)
}

rke.add <- function(rke1, rke2, verbose=F) {
  # Combines two RKE objects into one.
  # warning("rke.add() is not unit-tested")
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
  rke.all <- rke.update.new.pairs(rke=rke.all, keep.edges=all.edges)
  return(rke.all)
}

rke.keep.pairs = function(rke, pair.ids) {
  # Keeps the specified ids and removes all others.
  # Returns a <rke> object.
  all.pairs = rke.pair.ids(rke)
  rm.pairs = setdiff(all.pairs, pair.ids)
  return(rke.remove.pairs(rke, rm.pairs))
}

rke.subgraph <- function(rke, pair.type) {
  CHECK_MEMBER(pair.type, kPairTypes, msg="Correct pair type?")
  rke.keep.pairs(rke, pair.ids=rke.filter.pairs(rke, attr="pair.type", value=pair.type))
}

rke.edge.ids = function(rke) {
  if (rke.size(rke) == 0)
    return(c())
  as.vector(subset(rke$edges, can.donate==1)$edge.id)
}

rke.pair.ids = function(rke) { 
  as.vector(rke$pairs$pair.id)
}

rke.pairs.hospitals <- function(rke, pair.ids) {
  CHECK_UNIQUE(pair.ids)
  subrke = subset(rke$pairs, pair.id %in% pair.ids)
  return(subrke$hospital)
}

rke.hospital.pair.ids <- function(rke, hospital.id) {
  # Returns the pair ids that belong to that hospital
  return (rke$pairs$pair.id[which(rke$pairs$hospital == hospital.id)])
}

rke.cycles <- function(rke, include.3way=F) {
  # Computes 2-way (and) 3-way cycles for a specific rke object
  #
  # Returns:
  #   A data-frame with (type, id1, id2, id3) indicating that
  #   this is a "type"-way exchange id1->id2->id3
  #   If type=2 then id3=NA and if type=3 then it is assumed id3->id1 
  A2 = rke.2way.cycles(rke)
  A3 <- matrix(0, nrow=0, ncol=3)
  if (include.3way) 
    A3 <- matrix(rke.3way.cycles(rke), ncol=3)
  out = list(type=rep(0, 0),
             pair.id1=rep(0, 0), pair.id2=rep(0, 0), pair.id3=rep(0, 0))
  num.2cycles <- nrow(A2)
  num.3cycles <- nrow(A3)
  num.all <- num.2cycles + num.3cycles
  if (num.all > 0) {
    out$type <- rep(2, num.2cycles)
    out$pair.id1 <- A2[,1]
    out$pair.id2 <- A2[,2]
    out$pair.id3 <- rep(0, num.2cycles)
    out$type = c(out$type, rep(3, num.3cycles))
    out$pair.id1 <- c(out$pair.id1, A3[,1])
    out$pair.id2 <- c(out$pair.id2, A3[,2])
    out$pair.id3 <- c(out$pair.id3, A3[,3])
  }
  return (as.data.frame(out))
}

rke.cycles.membership <- function(rke, rke.cycles) {
  # Returns a (pairs x cycles) matrix that contains membership of pairs in cycles.
  # i.e. Aij = 1 if pair i belongs in cycle j
  # Row names *have* to be the same as the pair ids
  # If no cycles, return an empty matrix.
  all.pairs = rke.pair.ids(rke)
  # return a pairs x cycles matrix
  out = matrix(0, nrow=length(all.pairs), ncol=nrow(rke.cycles))
  if (length(all.pairs) == 0 | nrow(rke.cycles) == 0)
    return(out)
  subcycles = subset(rke.cycles, select=c(pair.id1, pair.id2, pair.id3))
  for (i in 1:length(all.pairs))
    out[i, ] <- as.numeric(rowSums(subcycles == all.pairs[i]))
  rownames(out) <- all.pairs
  return(out)
}

rke.filter.pairs <- function(rke, attrs, values) {
  # Filters the <pairs> object by matching the column "attr" to "value"
  #
  # Args:
  #   attrs : Array of RKE property (field)
  #   values : Array of the attrs values (should have the same length)
  # Returns:
  #   Array of pair ids.
  CHECK_MEMBER(attrs, names(rke$pairs))
  CHECK_TRUE(length(attrs) > 0, msg="Need some attributes")
  CHECK_EQ(length(attrs), length(values), msg="Same length of attr, value")
  if(rke.size(rke) == 0) {
    warning("Empty RKE")
    return(c())
  }
  matched.index <- 1:rke.size(rke)
  for (i in 1:length(attrs))
    matched.index <- intersect(matched.index, which(rke$pairs[[attrs[i]]] == values[i]))
  return(rke$pairs$pair.id[matched.index])
}

rke.2way.cycles <- function(rke) {
  # Returns a K x 2 matrix of 2-way cycles in the rke object.
  # (K = #2-way cycles)
  CHECK_rke(rke)
  out = matrix(0, nrow=0, ncol=2)
  if (rke.size(rke) == 0)
    return(out)
  A = rke$A  #  the adjacency matrix
  pair.ids = rownames(A)
  CHECK_EQ(pair.ids, rke.pair.ids(rke), msg="Pair ids from A and $pair.id should match")
  cycles.two = apply(which((A * t(A)) == 1, arr.ind=T), 2, function(i) as.numeric(pair.ids[i]))
  if (length(cycles.two) > 0) {
    x = apply(cycles.two, 1, function(r) r[2] > r[1])  # remove duplicates
    out = matrix(cycles.two[which(x), ], ncol=2)
  }
  return(out)
}

rke.3way.cycles <- function(rke) {
  # Given an adjacency matrix, it will return the 3-way cycles.
  #
  # Args: A adjacency matrix
  #
  # Returns: k x 3 matrix of cycles. Each row is (i,j,k) indicates i->j->k->i cycle
  CHECK_rke(rke)
  cycles.3way = matrix(0, nrow=0, ncol=3)
  if (rke.size(rke) == 0)
    return(cycles.3way)
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
  if (nrow(cycles.3way) == 0)
    return(cycles.3way)
  x = t(apply(cycles.3way, 1, diff))
  x = apply(x, 1, function(r) all(r < 0) || all(r > 0))
  ret = matrix(cycles.3way[which(x), ], ncol=3)
  # map back to pair ids.
  pair.ids = as.numeric(rownames(A))
  return (apply(ret, 2, function(i) pair.ids[i]))
}

rke.edge.by.pair <- function(rke, id.frame) {
  # Returns the list of edge id for the ids specified.
  #  Throws exception if not the correct format (e.g. unequal lengths)
  # 
  # Args:
  #   id.frame = data.frame with two columns (from->to)
  CHECK_EQ(class(id.frame), "data.frame", msg="Should be a data frame")
  CHECK_EQ(ncol(id.frame), 2, msg="from/to ids only")
  CHECK_MEMBER(0, id.frame, msg="No 0 pair ids")
  colnames(id.frame) <- c("pair.id1", "pair.id2")
  CHECK_DISJOINT(id.frame$pair.id1, id.frame$pair.id2, msg="No self-loops")
  subedges = subset(rke$edges, can.donate == 1)
  x1 = as.numeric(subedges$pair.id1 %in% id.frame$pair.id1)
  x2 = as.numeric(subedges$pair.id2 %in% id.frame$pair.id2)
  out = subedges$edge.id[which(x1 * x2 == 1)]
  CHECK_EQ(length(out), nrow(id.frame), msg="Not missing any edges")
  return (out)
}

plot.rke = function(rke, vertex.size=20) {
  library(igraph)
  if(rke.size(rke) == 0) {
    loginfo("Empty RKE can't be plotted.")
    return(TRUE)
  }
  g = graph.adjacency(rke$A, mode="directed")
  V(g)$color = rke$pairs$pair.color
  # igraph orders edges by source node.
  E(g)$color = arrange(subset(rke$edges, can.donate==1), pair.id1)$edge.color
  V(g)$label = str_c("H", rke$pairs$hospital, "#", rke$pairs$desc, rke$pairs$pair.id)
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,1))
  plot.igraph(g,layout=layout.auto, vertex.size=vertex.size)
}

# to-functions
pairid.to.pc <- function(rke, pair.ids) {
  if (rke.size(rke) == 0) return(c())
  CHECK_UNIQUE(pair.ids, msg="Pair ids need to be unique.")
  x = subset(rke$pairs, pair.id %in% pair.ids)
  return (x$pc)
}

pairid.to.hospital <- function(rke, pair.ids) {
  if (rke.size(rke) == 0) return(c())
  CHECK_UNIQUE(pair.ids, msg="Pair ids need to be unique.")
  x = subset(rke$pairs, pair.id %in% pair.ids)
  return (x$hospital)
}