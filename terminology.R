# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)
source("testing.R")
library(plyr)
library(logging)
library(stringr)
basicConfig()
# Terminology:
# (1) A "blood-type" is in {"O", "A", "B", "AB"} and represents a blood-type
#     A "blood-code" is an integer representation, in {1, 2, 3, 6} resp.
#     A "pair-type" is in {"O", "U", "S", "R"} and represent the over-demanded, 
#     under-demanded, self-demanded and reciprocal pairs respectively.
#     A "pair-code" is simply an integer from 1-16 which represents 
#     the 4x4 = 16 total blood-type dyads.
# (2) A "pairs" object is a DF that represents a collection of donor-patients
#     These are completely defined by (PairCodes, pra, hospital)
#     "pair codes" is a dataframe that describes the pairs themselves using 
#     (i) donor, patient types, (ii) blood.compatibility, (iii) pair as string
#     (iv) marginal blood-type probabilities. For example, see "kPairs"
#     that contains all such values for all 16 possible pair codes.
#  (3) "edges" is a DATAFRAME that contains id1, id2 vector that represent
#       pair dyads, and compatibility relationships from pair1->pair2.
#       pair.id1, pair.id2, blood.compatible, pra.compatible, can.donate
#         2          5          1                 1               1
#         5          8          1                 0               0
#       This means that pair 2 can donate to 5 but 5 cannot donate to 8 because
#       patient of 8 is PRA-sensitive to donor of 5.
# (4) An RKE ("random kidney exchange") defines a multi-hospital exchange pool.
#     It is a LIST of a <pairs> and <edges> objects. It also contains an adjacency
#     matrix Aij, where Aij = 1 if pair i can donate to pair j.
kAccuracy = 10^5
kBloodTypes  <- c("O", "A", "B", "AB")
kBloodCodes  <- c(1, 2, 3, 6)
kBloodTypeDistribution <- c(50, 30, 15, 5) / 100
get.blood.code.prob <- function(blood.code) {
  CHECK_MEMBER(blood.code, kBloodCodes, "Checking correct blood code")
  return(sapply(blood.code, function(i) kBloodTypeDistribution[which(kBloodCodes == i)]))
}
#   From type (e.g. "O") to numeric code
as.blood.code <- function(blood.type) {
  CHECK_MEMBER(blood.type, kBloodTypes)
  sapply(blood.type, function(bt) as.vector(kBloodCodes[which(kBloodTypes==bt)]))
}
as.blood.type <- function(blood.code) {
  CHECK_MEMBER(blood.code, kBloodCodes)
  sapply(blood.code,
         function(bc) as.vector(kBloodTypes[which(kBloodCodes==bc)]))
}
kUniformPRA <- 0.2
kNonUniformPRA <- c(0.05, 0.45, 0.9)
kNonUniformPRADistribution <- c(0.7, 0.2, 0.1)
kPairTypes <- c("R" ,"U", "O", "S")
kPairTypeColors <- list(R="yellow", U="gray", O="green", S="cyan")
# PAIRS object
# kPairs = 16 x 7 matrix:  Basic structure
# pc  donor  patient prob blood-type compatible     str  pair.type
# 1    1      1      0.25      1                   O-O     S
# 2    2      1      0.15      0                   A-O     U
#                 ...
kPairs <- expand.grid(donor=kBloodCodes, patient=kBloodCodes)
kPairs <- cbind(kPairs,
                     prob=apply(kPairs, 1, 
                                function(x) prod(kBloodTypeDistribution[which(kBloodCodes==x[1:2])])))
kPairs$prob <- as.vector(apply(
    apply(subset(kPairs, select=c(donor, patient)), 2, get.blood.code.prob),
    1, prod))
kPairs$blood.compatible <- as.numeric(with(kPairs, patient %% donor == 0))
kPairs$symmetric.compatible <- as.numeric(with(kPairs, donor %% patient == 0))
kPairs <- cbind(kPairs, desc=
                  apply(kPairs, 1, function(x) paste(as.blood.type(x[1:2]), collapse="-")))
kPairs$pair.type <- kPairTypes[1+ with(kPairs, 2 *blood.compatible + symmetric.compatible)]
kPairs$symmetric.compatible <- NULL
kPairs<- cbind(pc=1:16, kPairs)
kPairs$pair.color <- laply(kPairs$pair.type, function(i) kPairTypeColors[[i]])

rpra <- function(n, is.uniform=T) {
  ## Sample a PRA value. Uniform returns the same constant value (~0.2)
  if(is.uniform)
    return(rep(kUniformPRA, n))
  return(sample(kNonUniformPRA, size=n,
                replace=T, prob=kNonUniformPRADistribution))
}

rpairs <- function(n, pair.ids,
                   uniform.pra,
                   hospital.id=1,
                   blood.type.distr=kBloodTypeDistribution) {
  # Samples pairs for a kidney-exchange.
  # 
  # Args:
  #   n = #pairs to sample
  #   uniform.pra = T if PRA=constant or F if PRA is to be drawn from distribution.
  #   blood.type.distr = the blood-type distribution to use.
  # Returns:
  #   A <pairs> object. Recall this is "PairCodes" + PRAs  + hospital info 
  sample.pairs <- empty.pairs()
  sample.pras = c()
  while(nrow(sample.pairs) < n) {
    no.samples =  n + 100
    ## Sample the pair codes.
    pair.codes = sample(kPairs$pc, size=no.samples, replace=T, prob=kPairs$prob)
    ## Sample their PRAs
    pras =  rpra(no.samples, is.uniform=uniform.pra)
    ## Pair-internal PRA compatibilities.
    pra.compatible = rbinom(no.samples, size=1, prob=1-pras)
    ##  Pair-internal blood-type compatibility.
    blood.compatible = kPairs[pair.codes, ]$blood.compatible
    # enter the pool only if blood incompatible, or pra incompatible.
    enter.pool = which(blood.compatible * pra.compatible == 0)
    if(length(enter.pool) > 0) {
      sample.pairs = rbind(sample.pairs, kPairs[pair.codes[enter.pool], ])
      sample.pras = c(sample.pras, pras[enter.pool])
    }
  }
  CHECK_EQ(length(sample.pras), nrow(sample.pairs), "#Sampled PRAs=#pairs")
  subset.ids <- sample(nrow(sample.pairs), size=n, replace=F)
  sample.pairs = sample.pairs[subset.ids, ]
  rownames(sample.pairs) <- 1:nrow(sample.pairs)
  sample.pairs$hospital = c(sample.pairs$hospital, rep(hospital.id, n))
  sample.pairs$pra = sample.pras[subset.ids]
  sample.pairs$pair.id <- pair.ids
  return(sample.pairs)
}

# RKE definitions
empty.rke <- function() {
  obj = list(pairs=empty.pairs(),
             compatibility=data.frame(),
             A=matrix(0, nrow=0, ncol=0))
  class(obj) <- "rke"
  return(obj)
}
empty.pairs <- function() subset(kPairs, desc=="non.exists")
empty.edges <- function(size) {
  empty.data = rep(0, size)
  obj = data.frame(
    pair.id1=empty.data, 
    pair.id2=empty.data,
    blood.compatible=empty.data,
    pra.compatible=empty.data,
    self.loop=empty.data,
    can.donate=empty.data)
}

CHECK_rke <- function(rke) {
  CHECK_MEMBER(c("pairs", "edges", "A"), names(rke), "Match RKE fields.")
  CHECK_EQ(nrow(rke$pairs), nrow(rke$A))
  CHECK_MEMBER(as.vector(rke$A), c(0,1), msg="Aij in {0,1}")  # A is binary adjacency matrix.
  CHECK_SETEQ(diag(rke$A), c(0))  # no self-loops
  CHECK_INTERVAL(rke$pairs$prob, min=0, max=1, "Correct probabilities")
  CHECK_EQ(length(unique(rke$pairs$pair.id)), nrow(rke$A), "Pairs should have unique ids")
  CHECK_EQ(sum(rke$A), sum(rke$edges$can.donate),
           "Equal #edges in A and EDGES structs.")
  CHECK_pairs(rke$pairs)
  CHECK_edges(rke$edges)
}

CHECK_pairs <- function(pairs) {
  CHECK_MEMBER(c("pair.id", "donor", "patient", "pc", "hospital"), names(pairs))
  CHECK_MEMBER(pairs$pair.type, kPairTypes, "Correct pair types.")
  CHECK_TRUE(all(!duplicated(pairs$pair.id)), "No duplicate pair ids.")
}

CHECK_edges <- function(edges) {
  CHECK_MEMBER(c("pair.id1", "pair.id2", "can.donate"), y=names(edges))
  num.edges1 = length(unique(edges$pair.id1))
  CHECK_EQ(nrow(edges), num.edges1^2, msg="Cover all (pair1, pair2) combinations")
  CHECK_SETEQ(edges$pair.id1, edges$pair.id2, msg="Set(pairs1)=Set(pairs2)")
  CHECK_TRUE(all(!duplicated(edges$edge.id)), "No duplicate edge ids.")
}

sample.pairs.edges <- function(pairs, respect.edges=data.frame(), verbose=F) {
  # Given pairs, it will compute cross-compatilibities, based on blood-types
  # and cross-match PRA compatibilities.
  # Algorithm goes as follows:
  #   find all combinations of pairs (pair1, pair2)
  #   check whether donor1 -> patient2 (blood compatibility)
  #   check whether patient2 <- donor1 (PRA compatibility)
  num.pairs = nrow(pairs)
  # id1 = 1, 2, 4,..n, 1,2,3,...n, ...
  # id2 = 1,1,1,1,..,2,2,2,2....
  edges.row = data.frame(id1=rep(c(1:num.pairs), num.pairs),
                         id2=as.vector(sapply(1:num.pairs, function(i) rep(i, num.pairs))))
  num.dyads <- nrow(edges.row)
  if (verbose) loginfo("Grid expanded.")
  edges = empty.edges(num.dyads)
  if (verbose) loginfo("Created empty edges")
  edges$pair.id1 = pairs$pair.id[edges.row$id1]
  edges$pair.id2 = pairs$pair.id[edges.row$id2]
  if (verbose) loginfo("Grid for id rows expanded.")
  tmp.d1 = pairs$donor[edges.row$id1]  # donor of pair 1
  tmp.p2 = pairs$patient[edges.row$id2]  # patient of pair 2
  tmp.pra2 = pairs$pra[edges.row$id2]  # PRA of patient of 2
  if (verbose) loginfo("PRA computed.")
  # whether donor1 -> patient2 in terms of blood-type
  edges$blood.compatible = as.numeric(tmp.p2 %% tmp.d1==0)
  if (verbose) loginfo("Blood-type compatibility computed.")
  # whether donor1 -> patient2 in terms of PRA
  edges$pra.compatible = rbinom(num.dyads, size=1, prob=(1-tmp.pra2))
  if (verbose) loginfo("PRA compatilibity computed.")
  edges$self.loop = as.numeric(with(edges, pair.id1==pair.id2))
  #as.numeric(edges$pair.id1 == edges$pair.id2)
  if (verbose) loginfo("Self-loop computed.")
  edges$edge.id = with(edges, as.integer(kAccuracy * (log(2) * pair.id1 + log(3) * pair.id2)))
  #apply(subset(edges, select=c(pair.id1, pair.id2)), 1, function(x) paste(x, collapse=":"))
  if (verbose) loginfo("Edge-id computed..")
  edges$can.donate = with(edges, blood.compatible * pra.compatible * (1-self.loop))
  #rownames(edges) <- 1:num.dyads
  if (verbose) loginfo("Donation variable computed")
  rm(list=c("tmp.d1", "tmp.p2", "tmp.pra2"))
  if (verbose) loginfo("Removed some fields")
  if (nrow(respect.edges) > 0) {
    # need to respect the defined edges.
    matched.edge.ids = intersect(edges$edge.id, respect.edges$edge.id)
    edges.match.rows <- which(edges$edge.id %in% matched.edge.ids)
    respect.match.rows <- which(respect.edges$edge.id %in% matched.edge.ids)
    if (verbose) loginfo("Matched rows computed.")
    CHECK_EQ(length(respect.match.rows), length(edges.match.rows),
                    msg="Find all respected edges")
    if (length(edges.match.rows) > 0)
      edges[edges.match.rows, ] <- respect.edges[respect.match.rows, ]
  }
  if (verbose) loginfo("Respected earlier edges")
  return(edges)
}

map.edges.A <- function(edges) {
  # Creates the adj matrix of this edges object.
  CHECK_edges(edges)
  num.pairs = length(unique(edges$pair.id1))
  A = matrix(edges$can.donate, nrow=num.pairs)
  CHECK_EQ(nrow(A), ncol(A), "Adj matrix A is square")
  CHECK_SETEQ(diag(A), c(0), msg("No self-loops"))
  rownames(A) = edges$pair.id1[1:num.pairs]
  colnames(A) = edges$pair.id1[1:num.pairs]
  return(A)
}

rrke <- function(n, pair.ids=1:n, hospital.id=1,
                 uniform.pra=T,
                 blood.type.distr=kBloodTypeDistribution,
                 verbose=F) {
  # Sample a RKE object
  #
  # Args: n = # pairs in the pool.
  #       hospital.id = which hospital they belong to
  #       uniform.pra = use uniform PRA or not
  #       blood.type.distr = distribution over blood types.
  # pair.ids will be unique for this hospital.
  CHECK_EQ(length(pair.ids), n, msg="Should have #pair.ids=n")
  rke = empty.rke()
  if(n==0) return(rke)
  CHECK_EQ(length(pair.ids), n, msg="#pair ids = n")
  # 1. Sample the pairss
  pairs = rpairs(n, hospital.id=hospital.id, uniform.pra=uniform.pra,
                     blood.type.distr=blood.type.distr, pair.ids=pair.ids)
  if (verbose) loginfo("Pairs were sampled.")
  edges = sample.pairs.edges(pairs, verbose=verbose)
  if (verbose) loginfo("Edges were sampled")
  # Create the RKE object.
  rke = empty.rke()
  rke$pairs = pairs
  # remove self-loops
  rke$edges = edges # subset(as.data.frame(edges), pair.id1 != pair.id2)
  # number edges. These will be edge ids.
  rke$A = map.edges.A(rke$edges)
  if (verbose) loginfo("Adj matrix was created")
  return(rke)
}

append.rke <- function(rke1, rke2, verbose=F) {
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
  rke.all$edges = sample.pairs.edges(rke.all$pairs, respect.edges=all.edges,
                                     verbose=verbose)
  rke.all$A = map.edges.A(rke.all$edges)
  return(rke.all)
}

mu.thm = function(n) {
  0.556 * n  -0.338 * sqrt(n)- 2
}