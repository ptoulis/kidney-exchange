This code can be used to perform research on Kidney Exchanges.
In this domain, a random graph is formed where each node is a patient/donor pair 
(i.e. a donor willing to donate to a patient but can't because they are incompatible), and each directed edge denotes compatibility between pairs. The overall goal is to create efficient exchanges and also make sure that participants (hospitals) have the incentive to share truthfully their donor/patient pairs.

More information on the model can be found in (Toulis and Parkes, 2013) -- see here http://www.people.fas.harvard.edu/~ptoulis/harvard-homepage/papers/pkd-toulis-parkes2013.pdf
This is also the paper where the code is introduced and extensively used (please cite that paper wrt to this provided code).

The "xCM" mechanisms can be found in this paper. The "Bonus" mechanism can be found in (Ashlagi and Roth, 2013) -- see reference in the above paper.

### INSTALLATION
In order to be able to run the code you will need R (http://www.r-project.org/)
Any version will be ok as long as it agrees with the Gurobi R library.
The steps will be the following

1. Download the current "kidney-exchange" source code.
2. Register for a free account in www.gurobi.com
3. Get a free Gurobi academic license (you will need a valid .edu email). For that 
   you should go to Downloads -> Licenses
4. Download the Gurobi software, again through Downloads -> Gurobi software and follow the installation instructions.
   That should be straightforward: simply unpack the archive and make sure to inform the environment variables 
   (e.g. GUROBI_HOME)
5. Next you will need to install the R Gurobi library. There are two options: (i) Use install.packages("gurobi") 
   directly on the R console (this didn't work for me, for the latest R version) or (ii) find the "R" folder in the
   main Gurobi folder and find a file like gurobi-<version>_R_x86..tar.gz.  You can manually install this by typing 
   on a console something like "R CMD INSTALL <filename_of_tar_gz>". I had to do this for the latest version of R.
6. Last, you will need to activate the license by using the "grbgetkey" command (this is also explained in the Gurobi 
   documentation)

### SYNOPSIS
The functionality provided, covers the following: sample donor/patient graphs, clear/match graphs, create pool of graphs, run a clearing mechanism on a pool, compute outcomes/statistics on matchings etc.

The following code, creates a random patient/donor graph with 20 pairs, plots the graph, and then computes a maximum matching; finally it reports how many matches (transplantations) were performed.
```
   rke <- rrke(20) 
   plot(rke)
   m <- max.matching(rke)
   no.matches <- get.matching.utility(m)
```

The following code creates and clears a "multi-hospital" kidney-exchanges market, allowing for 3-cycles i.e., 
pair1 -> pair2 -> pair1. It creates a pool of 4 hospitals, each sharing a graph with 20 patient/donor pairs,
and then lets one hospital (Hospital 3) "deviate" i.e., hide all pairs that can be matched internally.
It then runs "xCM" internally to clear the market and prints the results.
```
pool = rrke.pool(m=4, n=20, uniform.pra=T)  # uniform PRA model (see Toulis & Parkes, 2013)
kpd = kpd.create(pool, strategy.str="ttc", include.3way=T) // include.3way allows for 3-cycles
out = Run.Mechanism(kpd, "xCM", include.3way=T)  // runs xCM on this pool
no.matches = get.matching.utility(out)  
print(out$information)  // prints detailed matching information (e.g. matching types, summaries etc)
```


### FILE ORGANIZATION
Currently, the project is composed by the following files

* File ```terminology.R``` contains all definitions for the terms used in the project. The main concepts/terms are    ```blood-type```, ```pairs```, ```edges```, ```matching```, ```rke``` (random kidney exchange), ```kpd``` (kidney    paired donation), ```strategy``` and ```mechanism```. Terminology also contains several generic/helper functions.
* File ```rke.R``` contains functions to manipulate a ```rke``` object. This structure represents a patient/donor     graph including information about the blood-types, graph (nodes and edges), pair types etc.
* File ```matching.R``` contains the interface to Gurobi, and computes maximum matchings using 2- or 3-cycles.
* File ```mechanisms.R``` contains the definition of a "market" (```kpd```) which is basically a set of ```rke```     objects in order to represent a multi-hospital setting. It implements random matching ("rCM"), extended matching    ("xCM", Toulis & Parkes, 2013) and "Bonus" (Ashlagi & Roth, 2013).
* File ```experiments.R``` contains *all* experiments in the paper (Toulis & Parkes, 2013)

### Contact
There are several more directions for this work and research in kidney exchanges in general.
One of them is the analysis and experimentation on more elaborate models, such as the "non-uniform PRA" model. 
This means modeling the sensitivity of patients.

Feel free to contact us if you have more questions or you need help.
My contact info is here http://www.people.fas.harvard.edu/~ptoulis/harvard-homepage/index.html

