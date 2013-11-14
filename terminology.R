# Copyright 2013 Panos Toulis, David C.Parkes
# Author: Panos Toulis(ptoulis@fas.harvard.edu)

library(plyr)
library(logging)
library(stringr)
basicConfig()
# Terminology:
# *  A "blood-type" is in {"O", "A", "B", "AB"} and represents a blood-type
#     A "blood-code" is an integer representation, in {1, 2, 3, 6} resp.
#     A "pair-type" is in {"O", "U", "S", "R"} and represent the over-demanded, 
#     under-demanded, self-demanded and reciprocal pairs respectively.
#     A "pair-code" is simply an integer from 1-16 which represents 
#     the 4x4 = 16 total blood-type dyads.
#  * A "pairs" object is a DF that represents a collection of donor-patients
#     These are completely defined by (PairCodes, pra, hospital)
#     "pair codes" is a dataframe that describes the pairs themselves using 
#     (i) donor, patient types, (ii) blood.compatibility, (iii) pair as string
#     (iv) marginal blood-type probabilities. For example, see "kPairs"
#     that contains all such values for all 16 possible pair codes.
#      A "pairs" data-frame has
#     pc  donor  patient  prob  blood.compatible  desc  pair.type pair.color  hospital pra  pair.id
#     3    3       1    0.075         0          B-O       U         gray       1      0.2     7
#  * "edges" is a DATAFRAME that contains id1, id2 vector that represent
#       pair dyads, and compatibility relationships from pair1->pair2.
#       pair.id1, pair.id2, blood.compatible, pra.compatible, self.loop   can.donate  edge.color  edge.id
#         2          5          1          tr       1               0          1           black     78232
#         5          8          1                 0               0          0           black     77723
#       This means that pair 2 can donate to 5 but 5 cannot donate to 8 because
#       patient of 8 is PRA-sensitive to donor of 5.
# *   A "matching" represents a transplant graph:
#     match = an <edges> data frame of the pairs involved
#     status = e.g. "OK" denoting the status of the matching (e.g. "INF")
#     utility = #pairs matched (as an integer)
# *  A strategy is a LIST(report, hide) that are vector of pair id's.
#    (report U hide) = {all pair ids} and (report ^ hide) = {}
# *  An RKE ("random kidney exchange") defines a multi-hospital exchange pool.
#     It is a LIST of a <pairs> and <edges> objects. It also contains an adjacency
#     matrix Aij, where Aij = 1 if pair i can donate to pair j.
# *  An RKE pool represent the pool created by several hospital RKE's.
#     It is defined as a list with "rke.list" having the list of the RKE's of
#     hospitals, and "rke.all" having the aggregate of the RKEs.
#      {  rke.list=>list(rke),  rke.all=>rke  }
# * A KPD object (Kidney-Paired donations) represents a collection of two
#     RKE pools: the "real" one that corresponds to what hospitals *have*
#     and the "reported" one that is the pool created from the pairs
#     that hospitals report. In symbols, strategy + RKE = KPD
# * A mechanism M is a function (rke.pool -> matching) 
#     The format of matching is the output of the max.matching function:
#     (see above)
#
# *  A "simulation setup" object (SimSetup) is a LIST that defines
#    nhospitals=#hospitals, sizes=ARRAY of hospital sizes,
#    uniform.pra, include.3way and nsims=#simulations
#   *  A "mechanism evaluation" object (mechEval) is a LIST of the following:
#      LIST(baseline= LIST(mech1=>{utility=c(..), matchingInfo}, mech2...),
#           deviation= LIST(mech1=...))
#     The basic concept is to compare between mechanisms.
kAccuracy = 10^5
kBloodTypes  <- c("O", "A", "B", "AB")
kBloodCodes  <- c(1, 2, 3, 6)
kPairCodes <- 1:16
kBloodTypeDistribution <- c(50, 30, 15, 5) / 100
# warning("Using special blood-type distribution")
# kBloodTypeDistribution <<- c(0, 50, 50, 0) / 100

get.blood.code.prob <- function(blood.code) {
  # Gets the marginal frequency of a specific blood code
  # e.g. get.blood.code.prob(1) = 0.5 i.e, "O" has 50% freq in the population
  CHECK_MEMBER(blood.code, kBloodCodes, "Checking correct blood code")
  return(sapply(blood.code, function(i) kBloodTypeDistribution[which(kBloodCodes == i)]))
}

as.blood.code <- function(blood.type) {
  # Get the blood code of a given blood-type
  # e.g. as.blood.code("O") = 1
  CHECK_MEMBER(blood.type, kBloodTypes)
  as.vector(sapply(blood.type, function(bt) as.vector(kBloodCodes[which(kBloodTypes==bt)])))
}

as.blood.type <- function(blood.code) {
  # Blood code -> blood type
  # e.g. as.blood.type(2) = "A"
  CHECK_MEMBER(blood.code, kBloodCodes)
  as.vector(sapply(blood.code,
         function(bc) as.vector(kBloodTypes[which(kBloodCodes==bc)])))
}

# Definition of kPairs
# This is the data frame with all possible combinations of patient-donor.
# Also includes information about compatibility, pair probability etc.
kUniformPRA <- 0.2
kNonUniformPRA <- c(0.05, 0.45, 0.9)
kNonUniformPRADistribution <- c(0.7, 0.2, 0.1)
kPairTypes <- c("R" ,"U", "O", "S")
kPairTypeColors <- list(R="yellow", U="gray", O="green", S="cyan")
# PAIRS object
# kPairs = 16 x 7 matrix:  Basic structure
# pc  donor  patient prob blood-type compatible     str  pair.type  pair.color
# 1    1      1      0.25      1        1           O-O     S
# 2    2      1      0.15      0        0           A-O     U
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
kPairs<- cbind(pc=kPairCodes, kPairs)
kPairs$pair.color <- laply(kPairs$pair.type, function(i) kPairTypeColors[[i]])

compatibility.mask = kUniformPRA * kPairs$blood.compatible + (1 - kPairs$blood.compatible)
kPairs$marginal.prob <- kPairs$prob * compatibility.mask
kPairs$marginal.prob <- kPairs$marginal.prob / sum(kPairs$marginal.prob)

# Functions on kPairs.
pc.to.desc <- function(pcs) {
  as.character(sapply(pcs, function(argpc) subset(kPairs, pc==argpc)$desc))
}

pc.to.pair.type <- function(pcs) {
  sapply(pcs, function(argpc) subset(kPairs, pc==argpc)$pair.type)
}

pc.reciprocal <- function(argpc) {
  x = subset(kPairs, pc == argpc)
  return(subset(kPairs, donor==x$patient & patient==x$donor)$pc)
}

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
  #   pair.ids = vector of pair ids (should be unique)
  #   uniform.pra = T if PRA=constant or F if PRA is to be drawn from distribution.
  #   blood.type.distr = the blood-type distribution to use.
  # Returns:
  #   A <pairs> object. Recall this will be:
  #   "PairCodes" + hospital + pra + pair.id
  if (n <= 0) stop("Require #npairs > 0")
  CHECK_EQ(n, length(pair.ids), msg="#pair.ids = #samples")
  CHECK_UNIQUE(pair.ids, msg="Pair ids should be unique.")
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
             edges=(empty.edges(0)),
             A=matrix(0, nrow=0, ncol=0))
  class(obj) <- "rke"
  return(obj)
}

empty.pairs <- function() {
  x = subset(kPairs, desc=="non.exists")
  x$hospital <- rep(0,0)
  x$pair.id <- rep(0, 0)
  x$pra <- rep(0, 0)
  return(x)
}

empty.edges <- function(size) {
  empty.data = rep(0, size)
  obj = data.frame(
    pair.id1=empty.data, 
    pair.id2=empty.data,
    blood.compatible=empty.data,
    pra.compatible=empty.data,
    self.loop=empty.data,
    can.donate=empty.data,
    edge.color=empty.data,
    edge.id=empty.data)
  return(obj)
}

CHECK_rke <- function(rke) {
  CHECK_SETEQ(c("pairs", "edges", "A"), names(rke), "Match RKE fields.")
  CHECK_EQ(nrow(rke$pairs), nrow(rke$A))
  CHECK_MEMBER(as.vector(rke$A), c(0,1), msg="Aij in {0,1}")  # A is binary adjacency matrix.
  if (nrow(rke$A) > 0)
    CHECK_SETEQ(diag(rke$A), c(0))  # no self-loops
  CHECK_EQ(rownames(rke$A), rke$pairs$pair.id, msg="rownames A == pair.id?")
  warning("Muted test in CHECK_rke().")
  # CHECK_EQ(sum(rke$A), sum(rke$edges$can.donate), "Equal #edges in A and EDGES structs.")
  CHECK_pairs(rke$pairs)
  CHECK_edges(rke$edges, rke$pairs)
}

CHECK_rke_subset <- function(rke.smaller, rke.bigger) {
  # Tests whether the smaller RKE is a subset of the bigger one
  #
  loginfo("Checking rke SUBSET")
  warning("SUBSET RKE test generates 2500 checks at most. This might not be optimal")
  nMaxTests <- 500
  CHECK_rke(rke.smaller)
  CHECK_rke(rke.bigger)
  # Pairs <= pairs(bigger) -- subset
  CHECK_MEMBER(rke.pair.ids(rke.smaller), rke.pair.ids(rke.bigger))
  CHECK_TRUE(all(table(rke.smaller$pair.type) <= table(rke.smaller$pair.type)))
  # Check the edges
  
  get.Aij <- function(rke, id1, id2) {
    i = which(rownames(rke$A)==id1)
    j = which(rownames(rke$A)==id2)
    return(rke$A[i,j])
  }
  new.pairs = rke.pair.ids(rke.smaller)
  ntrials = min(nMaxTests, nrow(rke.bigger$edges))
  edge.ids <- sample(1:nrow(rke.bigger$edges), size=ntrials, replace=F)
  loginfo(sprintf("Checking edges compatibility (%d tests)", ntrials))
  pb = txtProgressBar(style=3, min=0, max=1)
  tmp <- c()
  for (i in 1:ntrials) {
    edge.row.id = edge.ids[i]
    tmp <- c(tmp, edge.row.id)
    # pick a random edge from the original graph
    p1 <- rke.bigger$edges[edge.row.id,]$pair.id1
    p2 <- rke.bigger$edges[edge.row.id,]$pair.id2
    new.edge = subset(rke.smaller$edges, pair.id1==p1 & pair.id2==p2)
    if (p1 %in% new.pairs & p2 %in% new.pairs) {
      # p1 and p2 pairs are in the new graph
      old.edge = rke.bigger$edges[edge.row.id,]$can.donate
      new.edge = new.edge$can.donate
      # Check that the edges are the same in both old/new RKE
      CHECK_EQ(new.edge, old.edge)
      # Check that the A matrix has the same info
      CHECK_EQ(get.Aij(rke.smaller, p1, p2), new.edge)
    } else {
      CHECK_EQ(nrow(new.edge), 0)
    }
    setTxtProgressBar(pb, value= i/length(edge.ids))
  }
  CHECK_UNIQUE(tmp)  # simple check that the sampled edges are different
}

CHECK_rke.list <- function(rke.list) {
  CHECK_TRUE(length(rke.list) > 0, msg="Rke.list should not be empty.")
  # 1. Check every RKE
  laply(rke.list, CHECK_rke)
  hids <- rke.list.hospital.ids(rke.list)  # this is checking for valid hospital ids.
}

CHECK_pairs <- function(pairs) {
  CHECK_MEMBER(c("pair.id", "donor", "patient", "pc", "pra", "hospital"), names(pairs))
  if(nrow(pairs) > 0) {
    CHECK_UNIQUE(pairs$pair.id)
    CHECK_MEMBER(pairs$pc, kPairCodes, msg="pc in 1:16")
    CHECK_MEMBER(pairs$donor, kBloodCodes)
    CHECK_MEMBER(pairs$patient, kBloodCodes)
    CHECK_INTERVAL(pairs$prob, min=0, max=1, "prob in [0,1]")
    CHECK_MEMBER(pairs$pair.type, kPairTypes, "Correct pair types.")
    CHECK_TRUE(all(pairs$pair.id > 0), msg="All pair ids should be > 0")
    CHECK_INTERVAL(pairs$pra, min=0, max=1)
    CHECK_MEMBER(pairs$hospital, seq(min(pairs$hospital), max(pairs$hospital)))
  }
}

CHECK_edges <- function(edges, pairs) {
  # pairs = those used to create the edges.
  CHECK_MEMBER(c("pair.id1", "pair.id2", "can.donate"), y=names(edges))
  CHECK_UNIQUE(edges$edge.id)
  CHECK_MEMBER(c(edges$pair.id1, edges$pair.id2), pairs$pair.id)
}

CHECK_rke.pool <- function(rke.pool) {
  CHECK_SETEQ(names(rke.pool), c("rke.list", "rke.all"))
  CHECK_rke(rke.pool$rke.all)
  CHECK_rke.list(rke.pool$rke.list)
}

CHECK_kpd <- function(kpd) {
  CHECK_SETEQ(names(kpd), c("reported.pool", "real.pool"))
  CHECK_rke.pool(kpd$reported.pool)
  CHECK_rke.pool(kpd$real.pool)
}

CHECK_strategy <- function(strategy, all.pair.ids) {
  if (length(names(strategy)) == 1) {
    CHECK_SETEQ(all.pair.ids, strategy[[names(strategy)[1]]], msg="cover all pairs")
  } else {
    CHECK_SETEQ(names(strategy), c("hide", "report"))
    CHECK_DISJOINT(strategy$hide, strategy$report, msg="Report^Hide=0")
    CHECK_SETEQ(union(strategy$hid, strategy$report), all.pair.ids, msg="Union=all")
  }
}

simple.simSetup <- function() {
  return(list(nhospitals=4,
              sizes=c(10),
              uniform.pra=T,
              include.3way=T,
              nsims=10))
}

CHECK_simSetup <- function(simSetup) {
  CHECK_MEMBER(names(simSetup), c("nhospitals", "sizes", "uniform.pra", "include.3way", "nsims"))
  CHECK_TRUE(is.logical(simSetup$uniform.pra))
  CHECK_TRUE(is.logical(simSetup$include.3way))
  CHECK_GT(simSetup$nhospitals, 0)
  CHECK_TRUE(is.numeric(simSetup$sizes, 0))
  CHECK_GT(length(simSetup$sizes), 0)
  CHECK_GT(simSetup$nsims, 0)
}


rke.list.hospital.ids <- function(rke.list) {
  hids = 1:length(rke.list)
  lapply(rke.list, function(rke) CHECK_MEMBER(rke$pairs$hospital, hids))
  return (hids)
}

generate.pairs.edges <- function(pairs, keep.edges=c(), verbose=F) {
  # Given pairs, it will compute cross-compatilibities, based on blood-types
  # and cross-match PRA compatibilities.
  # Algorithm goes as follows:
  #   find all combinations of pairs (pair1, pair2)
  #   check whether donor1 -> patient2 (blood compatibility)
  #   check whether donor1 -> patient2 (PRA compatibility)
  # This function will sample edges according to blood-type and PRA
  # but will also respect the "keep.edges" (<edges> object)
  # i.e. if i->j in keep edges then i->j also in the new edges.
  #
  # Args: pairs = a "pairs" object (see above)
  #       eges = an "edges" object to be kept fixed
  #
  # Returns:
  #   An "edges" object (see terminology)
  num.pairs = nrow(pairs)
  CHECK_pairs(pairs)
  if(num.pairs == 0)
    return(empty.edges(0))
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
  edges$edge.color <- rep("black", nrow(edges))
  #rownames(edges) <- 1:num.dyads
  if (verbose) loginfo("Donation variable computed")
  rm(list=c("tmp.d1", "tmp.p2", "tmp.pra2"))
  if (verbose) loginfo("Removed some fields")
  if (length(keep.edges) > 0) {
    # need to respect the defined edges.
    matched.edge.ids = intersect(edges$edge.id, keep.edges$edge.id)
    edges.match.rows <- which(edges$edge.id %in% matched.edge.ids)
    respect.match.rows <- which(keep.edges$edge.id %in% matched.edge.ids)
    if (verbose) loginfo("Matched rows computed.")
    CHECK_EQ(length(respect.match.rows), length(edges.match.rows),
                    msg="Find all respected edges")
    if (length(edges.match.rows) > 0)
      edges[edges.match.rows, ] <- keep.edges[respect.match.rows, ]
    if (verbose) loginfo("Respected earlier edges")
  }
  # remove non-edges
  return(edges)
}

map.edges.adjacency <- function(edges, all.pair.ids) {
  # Creates the adj matrix of this edges object.
  # This is a n x n matrix, where n = #all.pair.ids
  CHECK_UNIQUE(all.pair.ids)
  n = length(all.pair.ids)
  A <<- matrix(0, nrow=n, ncol=n)
  subset.edges = subset(edges, can.donate == 1)
  get.pairIdIndex <- function(id) which(all.pair.ids == id)
  ddply(subset.edges, .(pair.id2),
        function(inlinks) {
          to.id = get.pairIdIndex(inlinks$pair.id2[1])
          in.ids = sapply(inlinks$pair.id1, get.pairIdIndex)
          A[in.ids, to.id] <<- 1
        })
  if (nrow(A) > 0)
    CHECK_SETEQ(diag(A), c(0), msg="No self-loops")
  CHECK_EQ(sum(A), sum(edges$can.donate), "A has correct #edges.")
  rownames(A) = all.pair.ids
  colnames(A) = all.pair.ids
  return(A)
}

rke.update.new.pairs <- function(rke, keep.edges) {
  # It will generate edges for this RKE object but will have
  # the "keep.edges" fixed. Should be called when "pairs" in the RKE object
  # are changed (e.g. when removing pairs, see remove.pairs in rke.R)
  edges = generate.pairs.edges(rke$pairs, keep.edges=keep.edges)
  rke$edges = edges  #subset(edges, can.donate == 1)
  # number edges. These will be edge ids.
  rke$A = map.edges.adjacency(rke$edges, rke$pairs$pair.id)
  return (rke)
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
  if (!is.numeric(pair.ids))
    stop("Pair ids need to be numeric for efficiency.")
  CHECK_EQ(length(pair.ids), n, msg="Should have #pair.ids=n")
  CHECK_UNIQUE(pair.ids, msg="Pair ids should be unique")
  rke = empty.rke()
  if(n==0) return(rke)
  # 1. Sample the pairss
  pairs = rpairs(n, hospital.id=hospital.id, uniform.pra=uniform.pra,
                     blood.type.distr=blood.type.distr, pair.ids=pair.ids)
  # Create the RKE object.
  rke = empty.rke()
  rke$pairs = pairs
  rke <- rke.update.new.pairs(rke, keep.edges=NULL)
  return(rke)
}

mu.thm = function(n, m=1) {
  0.556 * n *m -0.338 * sqrt(n * m)- 2
}

logthis <- function(x, verbose) {
  if (is.array(x))
    x = paste(x, collapse=", ")
  if (verbose) loginfo(x);
}

uniform.sample <- function(x) {
  if (length(x) == 0)
    stop("Non-empty vector needed for uniform sample.")
  if (length(x) == 1)
    return(x)
  sample(x, size=1)
}