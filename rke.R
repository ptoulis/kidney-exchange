# Panos Toulis, David C.Parkes
# 2012, Random Graph models for Kidney Exchanges
source("terminology.R")

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
 # "almost-equality" Check that the pair types have the same frequency table.
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
  remove.edges = subset(rke$edges, pair.id1 %in% rm.pair.ids | pair.id2 %in% rm.pair.ids)$edge.id
  rke$edges =  subset(rke$edges, !is.element(edge.id, remove.edges))
  rke$A = map.edges.adjacency(rke$edges, rke.pair.ids(rke))
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
  if(length(pair.ids) == 0)
    return(empty.rke())
  CHECK_TRUE(is.vector(pair.ids))
  CHECK_MEMBER(pair.ids, rke$pairs$pair.id)
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
  #
  #  Cycles = (type, id1, id2, id3)
  #             2     10   11   0
  #             3     20   24   50
  #                       ... 
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

rke.count.virtual.pairs <- function(xrke, return.bipartite=F) {
  # counts the number of virtual AB and BA pairs of the extended R subgraph.
  xrkeVirtualAB <- rke.keep.pairs(xrke, subset(xrke$pairs, desc %in% c("O-B", "A-O"))$pair.id)
  xrkeVirtualBA <- rke.keep.pairs(xrke, subset(xrke$pairs, desc %in% c("O-A", "B-O"))$pair.id)
  
  modify.xrke <- function(xrkeXY, Xdesc) {
    # This makes the extended-R-subgraph to be symmetric 
    # i.e. if OD->UD then UD->OD
    # We need this to decide how many virtual pairs are in this graph
    # through maximum matching.
    #
    # should be the over-demanded pair ids.
    #
    # check if no edges.
    if(sum(xrkeXY$edges$can.donate)==0)
      return(xrkeXY)
    pair.ids = subset(xrkeXY$pairs, desc==Xdesc)$pair.id
    # Xdesc = over-demanded pair.
    edgesXY = subset(xrkeXY$edges, can.donate==1 & pair.id1 %in% pair.ids)
    newEdgesXY = edgesXY
    # 1. Sample new edge ids for the symmetric edges.
    newEdgesXY$edge.id <- max(edgesXY$edge.id) + 100 + 1:nrow(newEdgesXY)
    # 2. Make edges symmetric
    newEdgesXY$pair.id1 <- edgesXY$pair.id2
    newEdgesXY$pair.id2 <- edgesXY$pair.id1
    xrkeXY$edges = rbind(edgesXY, newEdgesXY)
    xrkeXY$A = map.edges.adjacency(xrkeXY$edges, rke.pair.ids(xrkeXY))
    return(xrkeXY)
  }
  
  newAB = modify.xrke(xrkeVirtualAB, "O-B")
  newBA = modify.xrke(xrkeVirtualBA, "O-A")
  if(return.bipartite) {
    return(list(AB=newAB$pairs$pair.id, BA=newBA$pairs$pair.id))
  }
  nBA = get.matching.utility(max.matching(newBA)) / 2
  nAB = get.matching.utility(max.matching(newAB)) / 2
  return(list(AB=nAB, BA=nBA))
}

rke.extended.Rsubgraph <- function(rke) {
  # Returns the extended R subgraph
  # This includes R pairs and "virtual" pairs.
  #
  # TODO(ptoulis): Have a unit-test for this.
  pair.ids = subset(rke$pairs, desc %in% c("A-B", "B-A", "O-B", "B-O", "O-A", "A-O"))$pair.id
  xRke <- rke.keep.pairs(rke, pair.ids=pair.ids)
  xRke$edges <- subset(xRke$edges, can.donate==1)
  if(nrow(xRke$edges)==0) return(xRke)
  m1 = match(xRke$edges$pair.id1, xRke$pairs$pair.id)
  m2 = match(xRke$edges$pair.id2, xRke$pairs$pair.id)
  desc1 = xRke$pairs$desc[m1]
  H1 = xRke$pairs$hospital[m1]
  desc2 = xRke$pairs$desc[m2]
  H2 = xRke$pairs$hospital[m2]
  
  rm(rke)
  
  same.hospital = (H1==H2)
  
  vAB = (desc1=="O-B") & (desc2=="A-O") & same.hospital
  vBA = (desc1=="O-A") & (desc2=="B-O") & same.hospital
  u1 =  (desc1=="B-O")  &  (desc2=="A-B")
  u2 = (desc1=="A-O") & (desc2=="B-A")
  ab = (desc1=="A-B") & (desc2=="B-A" | desc2=="O-A")
  ba = (desc1=="B-A") & (desc2=="A-B" | desc2=="O-B")
  
  legit.edges = vAB | vBA | u1 | u2 | ab | ba
  
  xRke$edges <- xRke$edges[which(legit.edges), ]
  xRke$A = map.edges.adjacency(xRke$edges, rke.pair.ids(xRke))
  return(xRke)
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

plot.rke = function(rke, layout=layout.auto, 
                    mark.ex=10,
                    vertex.size=20,
                    font.size=1) {
  library(igraph)
  if(rke.size(rke) == 0) {
    loginfo("Empty RKE can't be plotted.")
    return(TRUE)
  }
  g = graph.adjacency(rke$A, mode="directed")
  V(g)$color = rke$pairs$pair.color
  V(g)$label.cex = font.size
  
  reverse.pair.names <- function(names) {
    out = c()
    for(name in names) {
      u = strsplit(name, split="")[[1]]
      dash.index = which(u=="-")
      i2 = seq(dash.index+1, length(u))
      i1 = seq(1, dash.index-1)
      new.name = str_c(paste(u[i2], collapse=""), "-",
                       paste(u[i1], collapse=""))
      out <- c(out, new.name)
    }
    return(out)
  }
  # igraph orders edges by source node.
  E(g)$color = arrange(subset(rke$edges, can.donate==1), pair.id1)$edge.color
  V(g)$label = str_c("H", rke$pairs$hospital, "#",
                     reverse.pair.names(as.character(rke$pairs$desc)),
                     rke$pairs$pair.id)
  
  par(mar=c(0,0,0,0))
  par(mfrow=c(1,1))
  groups = list()
  cols = c()

  for(group in unique(as.character(kPairs$pair.type))){
    groups[[group]] = subset(rke$pairs, pair.type==group)$pair.id
    colGroup = ifelse(group=="U", "white", subset(kPairs, pair.type==group)$pair.color[1])
    torgb = col2rgb(colGroup)[,1]/255
    cols <- c(cols, rgb(torgb[1], torgb[2], torgb[3], alpha=0.5))
  }
    print(groups)
  print(cols)
  plot.igraph(g, mark.groups=groups, mark.col=cols, 
              mark.border=cols, mark.expand=mark.ex,
              layout=layout, vertex.size=vertex.size)
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