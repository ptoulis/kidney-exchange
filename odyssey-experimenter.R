# Figure out current working directory
source("experiments.R")
#args <- as.integer(commandArgs(trailingOnly = T))

###     Changing the current working directory
##  HOW-TO run:
##  odyssey-experimenter    [sims]
##  All different sizes to try.. Can get up to 1000 but let's try this one first.
ns <- c(50,100, 150,200, 300)#, 500)

experiments <- c("E2","E3","E4")
scenarios  <- c("A","B","C","D")
k.hospitals <- c(3,3,3) ## only three can be


## Save results
save.e2 <- function(n, res) {

    random.id =as.integer(runif(1) * 10^6)
    filename = sprintf("./experiments/E2/%d/res_%d.Rdata", n, random.id)
    dir.create( sprintf("experiments/E2/%d", n), recursive=T, showWarnings = FALSE)
    save(res, file=filename)
}
save.e3 <- function(scenario, n, res) {

    random.id =as.integer(runif(1) * 10^6)
    filename = sprintf("./experiments/E3/%s/%d/res_%d.Rdata", scenario, n, random.id)
    dir.create( sprintf("experiments/E3/%s/%d", scenario, n), 
                    recursive=T, showWarnings = FALSE)
    save(res, file=filename)
}
save.e4 <- function(scenario, n, res) {

    random.id =as.integer(runif(1) * 10^6)
    filename = sprintf("./experiments/E4/%s/%d/res_%d.Rdata", scenario, n, random.id)
    dir.create( sprintf("experiments/E4/%s/%d", scenario, n), 
                recursive=T, showWarnings = FALSE)
    save(res, file=filename)
}


from.scenario.to.deviating <- function(k, scen) {
    if(scen=="A") return(c())
    if(scen=="B") return(c(1))
    if(scen=="C") return(c(2:k))
    if(scen=="D") return(c(1:k))
    stop("Not valid scenario")
}

Run.Sims = function(sims=10) {
    ##  Run Simulations
    for(i in 1:sims) {
        
        exp = sample(experiments,1)
        scenario = sample(scenarios, 1)
        k = sample(k.hospitals,1)
        n = sample(ns, 1)
        print(sprintf("[%d/%d] Running experiment %s k=%d n=%d scenario=%s",i, sims,exp, k, n, scenario))
        
        if(exp=="E2") {
            res = e2.once(n)
            random.id =as.integer(runif(1) * 10^6)
            filename = sprintf("./samples/%d/rke_%d.Rdata", n, random.id)
            dir.create( sprintf("samples/%d", n), recursive=T, showWarnings = FALSE)
            
            save.e2(n, res)
            
        }   else if(exp=="E3") {
            mech="rCM"
            res = e34.one.run(k, n, 
                              from.scenario.to.deviating(k, scenario), 
                              mech)
            
            save.e3(scenario, n, res)
            
        }  else if(exp=="E4") {
            mech="xCM"
            res = e34.one.run(k, n, 
                              from.scenario.to.deviating(k,scenario), 
                              mech)
            
            save.e4(scenario, n, res)
        }
    }
}

Run.Tables = function(ns.out, sims) {
    print("Running Experiment 1")
    out1 = E1(ns.out=ns.out, sims=sims)
    e1.to.latex(y=out1)    
    
    print("Running Experiment 2")
    out2 = E2(ns.out, sims)
    e2.to.latex(out2)
    
    run.exp = function(id, sc) {
        print(sprintf("Running Experiment %d -- Scenario %s", id, sc))
        out= NA
        if(id==3)
            out = E3(k=3, ns.out, sims, scenario=sc)
        else
            out = E4(k=3, ns.out, sims, scenario=sc)
        
        e34.to.latex(y=out, tableDesc= sprintf("E%d - Scenario %s", id, sc))
    }
    
    scenarios =  c("A","B","C","D")
    for(id in c(3,4))
        for(scen in scenarios)
            run.exp(id, scen)
}

