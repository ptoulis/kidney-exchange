# Synopsis has code fragments that show how to use the software
# and how this is structured.
source("mechanisms.R")
rke.pool = rrke.pool(m=3, n=30, uniform.pra=T)
kpd = kpd.create(rke.pool, strategy.str="ttt", verbose=T)
Run.Mechanism(kpd, mech="rCM", include.3way=T, verbose=T)