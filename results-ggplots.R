##  To-be deleted.
OUT.DIR="/home/ptoulis/Dropbox/A/bench/Kidney Exchanges (GEB)/images/new"

get.png.filename = function(png.name) {
  sprintf("%s/%s", OUT.DIR, png.name)
}
make.ggplot = function(title, filename, labels) {
  load(filename)
  # assumes results.
  y = list(scenario=c(), size=c(), rel.gain=c())
  for (scen in 1:length(names(results))) {
    scen.name = names(results)[scen]
    result = results[[scen.name]]
    for (size in names(result)) {
      A = result[[size]]
      rel.gains = A[1,]/A[2,]
      for(rel.gain in rel.gains)
      {
        y$scenario = c(y$scenario, scen.name)
        y$size = c(y$size, as.numeric(size))
        y$rel.gain = c(y$rel.gain, rel.gain)
      }
    }
  }
  
  y = data.frame(y, cutoff=rep(1, length(y$rel.gain)))
  
  a <- ggplot(y, main="YES",aes(x=factor(size), y=rel.gain))
  a <- a+facet_grid(scenario ~ size)
  a <- a+ geom_boxplot(aes(fill=factor(scenario))) 
  a <- a+ ylim(0.5, 2.0) + geom_hline(yintercept=1.0, lty=2, col="red", lwd=0.6)
  a <- a + labs(title = title, y="Scenario", x="hospital size (#pairs)")
  a <- a+ theme(plot.title = element_text(size = rel(2)))
  a <- a + scale_fill_manual(name="Scenarios", 
                             values=c("green", "orange", "red", "cyan"),
                             labels=labels)
  return(a)
}


make.incentives.plots <- function() {
  labels = c("PRA, (c,t,t) vs (t,t,t)",
             "PRA, (c,c,c) vs (t,c,c)",
            "~PRA, (c,t,t) vs (t,t,t)",
             "~PRA, (c,c,c) vs (t,c,c)")
  filename = "experiments/experiments/mech-rCM-m3-trials200.Rdata"
  title = "rCM incentives"
  rcm = make.ggplot(title=title, filename=filename, labels=labels)
  ggsave(rcm, filename= get.png.filename("incentives-rcm.png"))
  
  labels = c("PRA, (c,t,t) vs (t,t,t)",
             "PRA, (r,t,t) vs (t,t,t)",
             "~PRA, (c,t,t) vs (t,t,t)",
             "~PRA, (r,t,t) vs (t,t,t)")
  filename = "experiments/experiments/mech-xCM-m3-trials200.Rdata"
  title = "xCM incentives"
  xcm = make.ggplot(title=title, filename=filename, labels=labels)
  ggsave(xcm, filename=get.png.filename("incentives-xcm.png"))
  
  filename = "experiments/experiments/mech-Bonus-m4-trials200.Rdata"
  title = "Bonus incentives"
  bonus = make.ggplot(title=title, filename=filename, labels=labels)
  ggsave(bonus, filename=get.png.filename("incentives-bonus.png"))

}

get.efficiency.plot <- function(title, filename, labels) {
  load(filename)
  # assumes results.
  y = list(mech=c(), size=c(), rel.gain=c())
  for (mech in names(results)) {
    result = results[[mech]]
    if (mech=="rCM_t")
      next
    for (size in names(result)) {
      A = result[[size]]
      rel.gains = A / results[["rCM_t"]][[size]]
      for(rel.gain in rel.gains)
      {
        y$mech = c(y$mech, mech)
        y$size = c(y$size, as.numeric(size))
        y$rel.gain = c(y$rel.gain, rel.gain)
      }
    }
  }

  y = data.frame(y, cutoff=rep(1, length(y$rel.gain)))
  a <- ggplot(y, main="YES",aes(x=factor(size), y=rel.gain))
  a <- a+facet_grid(mech ~ size)
  a <- a+ geom_boxplot(aes(fill=factor(mech))) 
  a <- a+ ylim(0.76, 1.0) + geom_hline(yintercept=1.0, lty=2, col="red", lwd=0.6)
  a <- a + labs(title = title, y="Mechanism", x="hospital size (#pairs)")
  a <- a+ theme(plot.title = element_text(size = rel(2)))
  a <- a + scale_fill_manual(name="Mechanisms", 
                             values=c("green", "orange", "red"),
                             labels=labels)
  return(a)
  
}

make.efficiency.ggplots <- function() {
  a = get.efficiency.plot(title="Efficiency of mechanisms, uniform PRA", 
                          filename="experiments/experiments/efficiency-pra-TRUE-m4-trials200.Rdata",
                          labels=c("Bonus with long-R deviation",
                                   "rCM with long-R deviation", 
                                   "xCM at truth"))
  ggsave(a, file=get.png.filename("efficiency-pra-true.png"))
  a = get.efficiency.plot(title="Efficiency of mechanisms, non-uniform PRA", 
                          filename="experiments/experiments/efficiency-pra-FALSE-m4-trials200.Rdata",
                          labels=c("Bonus with long-R deviation",
                                   "rCM with long-R deviation", 
                                   "xCM at truth"))
  ggsave(a, file=get.png.filename("efficiency-pra-false.png"))
}

make.efficiency.manyhospitals.ggplots = function() {
  files = list.files("experiments/experiments/", pattern="vary", full.names=T)
  y = list(mech=c(), m=c(), rel.gain=c())
  for (file in files) {
    load(file)
    current.m = NA
    m.match = regexec(text=file, pattern="m(\\d+)-")
    current.m = regmatches(file, m.match)[[1]][2]
    
    for(mech in names(results)) {
      if(mech=="rCM_t") next
      rel.gains = results[[mech]][["25"]] / results[["rCM_t"]][["25"]]
      for(rel.gain in rel.gains) {
        y$mech = c(y$mech, mech)
        y$rel.gain = c(y$rel.gain, rel.gain)
        y$m = c(y$m, current.m)
      }
    }
  }
  labels = c("Bonus with long-R deviation", "rCM with long-R attack", "xCM at truth")
  title = "Efficiency of mechanisms, varying #hospitals"
  y = data.frame(y, cutoff=rep(1, length(y$rel.gain)))

  y$m2 = factor(sort(as.numeric(levels(x$m)[x$m])))
  a <- ggplot(y, aes(x=m2, y=rel.gain))
  a <- a+facet_grid(mech ~ m2)
  a <- a+ geom_boxplot(aes(fill=factor(mech))) 
  a <- a+ ylim(0.7, 1.0) + geom_hline(yintercept=1.0, lty=2, col="red", lwd=0.6)
  a <- a + labs(title = title, y="Mechanism", x="hospital size (#pairs)")
  a <- a+ theme(plot.title = element_text(size = rel(2)))
  a <- a + scale_fill_manual(name="Mechanisms", 
                             values=c("green", "orange", "red"),
                             labels=labels)
  ggsave(a, file=get.png.filename("efficiency-pra-false-varyhospitals.png"))
}
