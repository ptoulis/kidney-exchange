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

rke.keep.pairs = function(rke, pair.ids) {
  all.pairs = rke.pair.ids(rke)
  rm.pairs = setdiff(all.pairs, pair.ids)
  return(rke.remove.pairs(rke, rm.pairs))
}

rke.edge.ids = function(rke) rke$edges$edge.id
rke.pair.ids = function(rke) rke$pairs$pair.id

plot.rke = function(rke, vertex.size=20) {
  library(igraph)
  g = graph.adjacency(rke$A, mode="directed")
  V(g)$color = rke$pairs$pair.color
  V(g)$label = str_c("H", rke$pairs$hospital, "#", rke$pairs$desc, rke$pairs$pair.id)
  par(mar=c(0,0,0,0))
  plot.igraph(g,layout=layout.auto, vertex.size=vertex.size)
}