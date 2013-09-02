# Synopsis has code fragments that show how to use the software
# and how this is structured.
rm(list=ls())
source("mechanisms.R")
Rprof(filename="debug/profile.txt", memory.profiling=T)
rke.pool = rrke.pool(m=3, n=30, uniform.pra=T)
kpd = kpd.create(rke.pool, strategy.str="ctt", verbose=T)
Run.Mechanism(kpd, mech="rCM", include.3way=T, verbose=T)
Rprof(NULL)