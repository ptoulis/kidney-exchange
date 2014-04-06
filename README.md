This code can be used to perform research on Kidney Exchanges.
In this domain, a random graph is formed where each node is a donor/patient
and each directed edge denotes compatibility between pairs.
The overall goal is to create efficient exchanges and also make sure that participants (hospitals) have the incentive to share truthfully their donor/patient pairs.

More information on the model can be found in (Toulis and Parkes, 2013) -- see here http://www.people.fas.harvard.edu/~ptoulis/harvard-homepage/papers/pkd-toulis-parkes2013.pdf

The "xCM" mechanisms can be found in this paper. The "Bonus" mechanism can be found in (Ashlagi and Roth, 2013) -- see reference in the above paper.

INSTALLATION
============
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


