-- README file ---


CODE
=====
"rke" =  basic object of a Random Kidney Exchange. Contains graph, list of pairs, edges etc.

rke.R  -- contains all functions to manipulate rke's
		  e.g. sample, max.matching, get.pairs.attribute, pool.rke(list)
		  SYNOPSIS:
		  ---------
		  rke = rrke(50)
		  m = max.matching(rke)
		  print(m$matching$matched.ids)
		  
		  ##  this should clog the machine.
		  rke.list = rrke.many(3, 1000)
		  rke.all = pool.rke(rke.list)

experiments.R -- contains code to run experiments: 
                 E2 = violations, E3 = rCM,  E4 = xCM
				 Scenarios :  A=all truthful (heaviest)  B = 1 deviates, C = only 1 truthful, D=all deviate
				 SYNOPSIS:
				 --------
				 E3(k=3, ns.out=10, sims=100, scenario="B) -- runs 3 hospitals, 20 to 200 pairs, 100 simulations, B scenario
				 
mechanisms.R -- contains code to run the mechanisms, 
				rCM and xCM



ODYSSEY
=======
odyssey-*.R  -- code to run on Odyssey.
The "Sampler" and the "Experimenter" are called as commands and so 
there are respective .bsub files in the "bsub" folder. 

"-sampler.R"      -- will sample RKE objects
                     Saves on the "samples" folder.
"-experimenter.R"  -- will perform a random experiment, on a predefined size.
                      Saves on the "experiments" folder. 
"-interface.R"     -- attempt to bridge with Python (under construction)


To work on Odyssey:
-------------------
1. Run Putty (load "Odyssey") and login. 
    Password starts with QA....

2. Start the "Odyssey-Auth" application (keylock on the desktop)
    Use this for the confirmation code.

3. When logged in ,run "ssh airoldi01"

4. Run the WinSCP for each file transfer



ISSUES:
=======
Very inefficient to 
1.  Sample big RKEs, e.g. rke = rrke(10000)

	## Current performance: Sampe rrke(300) => 1.53s
								  rrke(600) => 11s
	                              rrke(800) => 39s
								  rrke(1000) => 45s
								  
2.  Pool large rkes, e.g.   rke.list = rrke.many(3, 1000)
			                rke.all = pool.rke(rke.list)
	## example      rke.list = rrke.many(k=3, n=300)   => 4.69s
					pool =    pool.rke(rke.list)       = 31.54s
3.  Maximum matching, e.g.  m = max.matching(rke)
	Python code is also available in "./python/mwmatching.py" 
	but it is not complete since it cannot implement the "practical xCM", since it cannot consider IR constraints.




EXPERIMENTS
================
1.    Model assumptions, µ(n):  (Table 2.) -- (in tables.R)   table1(sizes=.., trials=..)    
2.    Violations:               (Table 3.) --  (in experiments.R)    E2(sizes=..., sims=...)    then e2.to.latex(..)









