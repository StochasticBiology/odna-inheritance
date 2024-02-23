library(ggplot2)
library(ggpubr)

logit = function(x) {
  xt = x
  return(log(xt/(1-xt)))
}

npop = 100

sf = 2

######### baseline behaviour -- no leakage or environmental change
scale = "0.500"; penalty = "0.000"; 
npop.comp.df = data.frame()

for(npop in c(50, 100, 200)) {
  
  # read in all parallelised outputs for a given scale/penalty combination
  df = data.frame()
  expt=0
  tdf = read.csv(paste0("inherit-mean-out-", npop, "-", expt, "-", scale, "-", penalty, ".csv"))
  df = rbind(df, tdf)
  
  df = df[!is.na(df$expt),]
  
  my.trans = function(lambda) {
    return( ifelse(lambda == 0, log2(0.002), log2(lambda)))
  }
  
  ##### first look at the long-time behaviour (**NOT guaranteed to be steady state!)
  long.df = df[df$DUI == 0 & df$t == 1995,]
  
  mean.df = long.df[long.df$expt==0,]
  mean.df$meanmean = 0
  mean.df$meanvar = 0
  for(i in 1:nrow(mean.df)) {
    mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                                long.df$nDNA == mean.df$nDNA[i] &
                                                long.df$mu == mean.df$mu[i] &
                                                long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
    mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
  }
  
  npop.comp.df = rbind(npop.comp.df, mean.df)
  
  g.baseline = ggplot(mean.df[mean.df$env==0,], aes(x=leakage,y=nDNA,fill=meanmean/nDNA)) + 
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    geom_tile() + facet_grid(~ mu)
  
  
  
  # empirical fit to fitness behaviour for constant environments
  fit.df = mean.df[mean.df$env == 0 & mean.df$meanmean/mean.df$nDNA != 1,]
  fit.df$x = 1- fit.df$mu*(fit.df$nDNA**0.5)/(1-fit.df$leakage)
  fit.df$y = fit.df$meanmean/fit.df$nDNA
  my.fit = summary(lm(y ~ x, data=fit.df))
  plot.title = paste0("Empirical fit for constant environments\nR^2 = ", round(my.fit$r.squared, digits=3))
  
  emp.fit.g = ggplot(fit.df, aes(x=x, y=y, color=factor(leakage))) + 
    geom_point() + #+ geom_line() 
    labs(x="1 - mu*sqrt(nDNA)/(1 - lambda)", y="meanf", color="lambda") +
    geom_abline() +
    ggtitle(plot.title)
  
  png(paste0("constant-env-and-fit-", npop, ".png"), width=800*sf, height=300*sf, res=72*sf)
  print (ggarrange(g.baseline, emp.fit.g, nrow = 1, widths=c(2,1), labels=c("A", "B")) )
  dev.off()
}

# awkward code for Npop = 50, 100, 200 comparison at baseline
npop.comp.2.df = npop.comp.df
npop.comp.2.df$npop.50 = npop.comp.2.df$npop.100 = npop.comp.2.df$npop.200 = NA
for(i in 1:nrow(mean.df)) {
  ref.50 = which(npop.comp.df$env == npop.comp.2.df$env[i] &
                   npop.comp.df$nDNA == npop.comp.2.df$nDNA[i] &
                   npop.comp.df$mu == npop.comp.2.df$mu[i] &
                   npop.comp.df$leakage == npop.comp.2.df$leakage[i] &
                   npop.comp.df$Npop == 50)
  ref.100 = which(npop.comp.df$env == npop.comp.2.df$env[i] &
                    npop.comp.df$nDNA == npop.comp.2.df$nDNA[i] &
                    npop.comp.df$mu == npop.comp.2.df$mu[i] &
                    npop.comp.df$leakage == npop.comp.2.df$leakage[i] &
                    npop.comp.df$Npop == 100)
  ref.200 = which(npop.comp.df$env == npop.comp.2.df$env[i] &
                    npop.comp.df$nDNA == npop.comp.2.df$nDNA[i] &
                    npop.comp.df$mu == npop.comp.2.df$mu[i] &
                    npop.comp.df$leakage == npop.comp.2.df$leakage[i] &
                    npop.comp.df$Npop == 200)
  if(length(ref.50) > 0 ) {  npop.comp.2.df$npop.50[i] = npop.comp.df$meanmean[ref.50]/npop.comp.df$nDNA[ref.50] }
  if(length(ref.100) > 0 ) { npop.comp.2.df$npop.100[i] = npop.comp.df$meanmean[ref.100]/npop.comp.df$nDNA[ref.100] }
  if(length(ref.200) > 0 ) { npop.comp.2.df$npop.200[i] = npop.comp.df$meanmean[ref.200]/npop.comp.df$nDNA[ref.200] }
  
}

ggplot(npop.comp.2.df) + geom_point(aes(x=npop.100, y=npop.50, color=factor(nDNA), shape="50"), size=3) +
  geom_point(aes(x=npop.100, y=npop.200, color=factor(nDNA), shape="200"), size=3) + geom_abline() +
  labs(x = "Mean fitness, Npop = 100", y = "Mean fitness, different Npop", color = "N", shape="Npop") +
  xlim(0.7, 1) + ylim(0.7, 1)

############# changing environments
###### phase space behaviour
# wrapper script currently only supports loop for npop = 100.
# for npop = 100, comment single vals and use loops; otherwise comment loops and use single values

#npop = 100
comps.df = data.frame()
scale = "0.500"; penalty = "0.000"; 
for(npop in c(50, 100, 200)) {{
  #  for(scale in c("0.000", "0.500")) {
  #    for(penalty in c("0.000", "1.000")) {
  
  # read in all parallelised outputs for a given scale/penalty combination
  df = data.frame()
  for(expt in 0:6) {
    tdf = read.csv(paste0("inherit-mean-out-", npop, "-", expt, "-", scale, "-", penalty, ".csv"))
    df = rbind(df, tdf)
  }
  df = df[!is.na(df$expt),]
  
  my.trans = function(lambda) {
    return( ifelse(lambda == 0, log2(0.002), log2(lambda)))
  }
  
  ##### first look at the long-time behaviour (**NOT guaranteed to be steady state!)
  long.df = df[df$DUI == 0 & df$t == 1995,]
  
  mean.df = long.df[long.df$expt==0,]
  mean.df$meanmean = 0
  mean.df$meanvar = 0
  for(i in 1:nrow(mean.df)) {
    mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                                long.df$nDNA == mean.df$nDNA[i] &
                                                long.df$mu == mean.df$mu[i] &
                                                long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
    mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
  }
  
  comps.df = rbind(comps.df, mean.df)
  
  g.dui0 = ggplot(mean.df, aes(x=leakage,y=nDNA,fill=meanmean/nDNA)) + 
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    geom_tile() + facet_grid(mu ~ env)
  
  g.baseline = ggplot(mean.df[mean.df$env==0,], aes(x=leakage,y=nDNA,fill=meanmean/nDNA)) + 
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    geom_tile() + facet_grid(~ mu)
  
  ##### first look at the long-time behaviour (**NOT guaranteed to be steady state!)
  long.df = df[df$DUI == 1 & df$t == 1995,]
  
  mean.df = long.df[long.df$expt==0,]
  mean.df$meanmean = 0
  mean.df$meanvar = 0
  for(i in 1:nrow(mean.df)) {
    mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                                long.df$nDNA == mean.df$nDNA[i] &
                                                long.df$mu == mean.df$mu[i] &
                                                long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
    mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
  }
  
  comps.df = rbind(comps.df, mean.df)
  
  g.dui1 = ggplot(mean.df, aes(x=leakage,y=nDNA,fill=meanmean/nDNA)) + 
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    geom_tile() + facet_grid(mu ~ env)
  
  print(paste0("Outputting ", scale, " ", penalty))
  
  sf = 2
  png(paste0("panels-", npop, "-", scale, "-", penalty, ".png"), width=800*sf, height=1000*sf, res=72*sf)
  print( ggarrange(g.dui0, g.dui1, nrow=2, ncol=1, labels = c("A. DUI = 0", "B. DUI = 1")) )
  dev.off()
  
  #### now look at the time series to see if steady state looks reasonable
  t.df = df[df$DUI == 0 & df$t < 1950,]
  g.traces = ggplot(t.df, aes(x=t, y=mean.f/nDNA, color=factor(leakage+nDNA))) + geom_point(size=0.5) +
    facet_grid(mu ~ env) + theme(legend.position = "none")
  
  sf = 2
  png(paste0("panels-all-", npop, "-", scale, "-", penalty, ".png"), width=1000*sf, height=1000*sf, res=72*sf)
  print( ggarrange( ggarrange(g.dui0+theme(legend.position="none"), g.dui1+theme(legend.position="none"), nrow=1, ncol=2, labels = c("A. DUI = 0", "B. DUI = 1")),
                    g.traces, nrow=2, ncol=1, labels=c("", "C. Time series") ) )
  dev.off()
  
}
}

# awkward code for comparisons
# mean.df currently holds the DUI=1 set of observations 
# we want the set of other fields from this, as is for the DUI comparison...
dui.comp.df = mean.df
dui.comp.df$DUI.1 = dui.comp.df$DUI.0 = NA
# ... and with two different DUI values for the Npop comparison
npop.comp.df = rbind(data.frame(mean.df, DUI=0), data.frame(mean.df, DUI=1))
npop.comp.df$npop.50 = npop.comp.df$npop.100 = npop.comp.df$npop.200 = NA

# go through this set of fields, pulling mean fitness for each param set with and without DUI
for(i in 1:nrow(dui.comp.df)) {
  ref.1 = which(comps.df$env == dui.comp.df$env[i] &
                  comps.df$nDNA == dui.comp.df$nDNA[i] &
                  comps.df$mu == dui.comp.df$mu[i] &
                  comps.df$leakage == dui.comp.df$leakage[i] &
                  comps.df$Npop == dui.comp.df$Npop[i] &
                  comps.df$DUI == 1)
  ref.0 = which(comps.df$env == dui.comp.df$env[i] &
                  comps.df$nDNA == dui.comp.df$nDNA[i] &
                  comps.df$mu == dui.comp.df$mu[i] &
                  comps.df$leakage == dui.comp.df$leakage[i] &
                  comps.df$Npop == dui.comp.df$Npop[i] &
                  comps.df$DUI == 0)
  
  dui.comp.df$DUI.1[i] = comps.df$meanmean[ref.1]/comps.df$nDNA[ref.1]
  dui.comp.df$DUI.0[i] = comps.df$meanmean[ref.0]/comps.df$nDNA[ref.0]
}

# go through this set of fields, pulling mean fitness for each param set with Npop = 50, 100, 200
for(i in 1:nrow(npop.comp.df)) {
  ref.50 = which(comps.df$env == npop.comp.df$env[i] &
                   comps.df$nDNA == npop.comp.df$nDNA[i] &
                   comps.df$mu == npop.comp.df$mu[i] &
                   comps.df$leakage == npop.comp.df$leakage[i] &
                   comps.df$DUI == npop.comp.df$DUI[i] &
                   comps.df$Npop == 50)
  ref.100 = which(comps.df$env == npop.comp.df$env[i] &
                    comps.df$nDNA == npop.comp.df$nDNA[i] &
                    comps.df$mu == npop.comp.df$mu[i] &
                    comps.df$leakage == npop.comp.df$leakage[i] &
                    comps.df$DUI == npop.comp.df$DUI[i] &                    
                    comps.df$Npop == 100)
  ref.200 = which(comps.df$env == npop.comp.df$env[i] &
                    comps.df$nDNA == npop.comp.df$nDNA[i] &
                    comps.df$mu == npop.comp.df$mu[i] &
                    comps.df$leakage == npop.comp.df$leakage[i] &
                    comps.df$DUI == npop.comp.df$DUI[i] &                 
                    comps.df$Npop == 200)
  if(length(ref.50) > 0 ) {  npop.comp.df$npop.50[i] = comps.df$meanmean[ref.50]/comps.df$nDNA[ref.50] }
  if(length(ref.100) > 0 ) { npop.comp.df$npop.100[i] = comps.df$meanmean[ref.100]/comps.df$nDNA[ref.100] }
  if(length(ref.200) > 0 ) { npop.comp.df$npop.200[i] = comps.df$meanmean[ref.200]/comps.df$nDNA[ref.200] }
  
}

dui.effect.g = ggplot(dui.comp.2.df, aes(x=DUI.0, y=DUI.1, color=factor(env))) + 
  geom_point(alpha=0.5) + geom_abline() + 
  labs(x="Mean fitness without DUI", y="Mean fitness with DUI", color="Environment")

npop.effect.g = ggplot(npop.comp.df) + geom_point(aes(x=npop.100, y=npop.50, color=factor(env), shape="50"), size=2, alpha=0.5) +
  geom_point(aes(x=npop.100, y=npop.200, color=factor(env), shape="200"), size=2, alpha=0.5) + geom_abline() +
  labs(x = "Mean fitness, Npop = 100", y = "Mean fitness, different Npop", color = "Environment", shape="Npop")

png("npop-dui-effect.png", width=700*sf, height=350*sf, res=72*sf)
print(ggarrange(npop.effect.g, dui.effect.g, labels=c("A. Npop influence", "B. DUI influence")))
dev.off()

########### zoomed-in region for simple case

scale = "0.500"
penalty = "0.000"
df = data.frame()
for(expt in c(0:10,20,50,100,200)) {
  tdf = read.csv(paste0("inherit-zoom-mean-out-", expt, "-", scale, "-", penalty, ".csv"))
  df = rbind(df, tdf)
}
df = df[!is.na(df$expt),]

df$period = as.numeric(df$env)

long.df = df[df$DUI == 0 & df$t > 1900,]

mean.df = long.df[long.df$expt==0,]
mean.df$meanmean = 0
mean.df$meanvar = 0
for(i in 1:nrow(mean.df)) {
  mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
  mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                            long.df$nDNA == mean.df$nDNA[i] &
                                            long.df$mu == mean.df$mu[i] &
                                            long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
}

g.zoom = ggplot(mean.df, aes(x=leakage,y=nDNA,fill=meanmean/nDNA)) + 
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
  geom_tile() + facet_wrap(~ period)

g.zoom

sf = 2
png(paste0("zoom-region.png"), width=800*sf, height=600*sf, res=72*sf)
print(g.zoom)
dev.off()


##### heteroplasmy penalty, different fitness cost

g.models = list()
npop = 100
for(model in 1:2) {
  if(model == 1) { scale = "0.250"; penalty = "0.000" }
  if(model == 2) { scale = "0.500"; penalty = "0.100" }
  
  # read in all parallelised outputs for a given scale/penalty combination
  df = data.frame()
  for(expt in 0:6) {
    tdf = read.csv(paste0("inherit-mean-out-", npop, "-", expt, "-", scale, "-", penalty, ".csv"))
    df = rbind(df, tdf)
  }
  df = df[!is.na(df$expt),]
  
  long.df = df[df$DUI == 0 & df$t == 1995,]
  long.df$mu[long.df$mu == 0] = 1e-6
  
  mean.df = long.df[long.df$expt==0,]
  mean.df$meanmean = 0
  mean.df$meanvar = 0
  for(i in 1:nrow(mean.df)) {
    mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                                long.df$nDNA == mean.df$nDNA[i] &
                                                long.df$mu == mean.df$mu[i] &
                                                long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
    mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
  }
  g.models[[length(g.models)+1]] = ggplot(mean.df, aes(x=leakage, y=nDNA, fill=meanmean/nDNA)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(mu ~ env)
  
}


ggarrange(plotlist = g.models)

#### redundant from here on except for confirming that f = 1-mu when N = 1
######## baseline behaviour for single oDNA per cell

scale = "0.500"; penalty = "0.000"
npop = 100; expt = 0

# read in all parallelised outputs for a given scale/penalty combination
df = data.frame()
tdf = read.csv(paste0("inherit-baseline-single-mean-out-", npop, "-", expt, "-", scale, "-", penalty, ".csv"))
df = rbind(df, tdf)
tdf = read.csv(paste0("inherit-baseline-mean-out-", npop, "-", expt, "-", scale, "-", penalty, ".csv"))
tdf$Npop=100
df = rbind(df, tdf)

df = df[!is.na(df$expt),]

long.df = df[df$DUI == 0 & df$t == 1995,]
long.df$mu[long.df$mu == 0] = 1e-6

mean.df = long.df[long.df$expt==0,]
mean.df$meanmean = 0
mean.df$meanvar = 0
for(i in 1:nrow(mean.df)) {
  mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
  mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                            long.df$nDNA == mean.df$nDNA[i] &
                                            long.df$mu == mean.df$mu[i] &
                                            long.df$leakage == mean.df$leakage[i]], na.rm= TRUE)
}


#log((1-mu)/(1-(1-mu)))
#= log((1-mu)/mu)
#= log(1/mu - 1)
# ~= log(1/mu)

# we have a choice of just matching the f > 0.5 behaviour -- which we csn do with a 1-f(mu, N) form 
# or matching the full sigmoid, which needs logit-like terms

### sigmoid picture
ggplot(mean.df[mean.df$meanmean/mean.df$nDNA != 1,]) + geom_point(aes(x=mu, y=logit(meanmean/nDNA), color=factor(nDNA))) +
  geom_line(aes(x=mu, y=-log(mu*nDNA**0.5), color=factor(nDNA))) +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ ylim(NA,5)

# transformed -- N=1 case not so well fit, but full sigmoid shape captured
ggplot(mean.df[mean.df$meanmean/mean.df$nDNA != 1,]) + geom_point(aes(x=mu, y=meanmean/nDNA, color=factor(nDNA))) +
  geom_line(aes(x=mu, y=1/(1+ mu*nDNA**0.5), color=factor(nDNA))) +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ ylim(NA,5)

# predicted vs observed
ggplot(mean.df[mean.df$meanmean/mean.df$nDNA != 1,]) + geom_point(aes(x=1/(1+ mu*nDNA**0.5), y=meanmean/nDNA, color=factor(nDNA))) +
  # geom_line(aes(x=mu, y=1/(1+ mu*nDNA**0.5), color=factor(nDNA))) +
  # scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ ylim(NA,5)

### 1-mu picture for f > 0.05  -- N=1 case better fit, but tails no good
ggplot(mean.df[mean.df$meanmean/mean.df$nDNA != 1,]) + geom_point(aes(x=mu, y=meanmean/nDNA, color=factor(nDNA))) +
  geom_line(aes(x=mu, y=1-mu*nDNA**0.5, color=factor(nDNA))) +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0.,1)

# predicted vs observed
ggplot(mean.df[mean.df$meanmean/mean.df$nDNA != 1,]) + geom_point(aes(x=1- mu*nDNA**0.5, y=meanmean/nDNA, color=factor(nDNA))) +
  # geom_line(aes(x=mu, y=1/(1+ mu*nDNA**0.5), color=factor(nDNA))) +
  # scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+ ylim(NA,5)

to.fit = mean.df[mean.df$meanmean/mean.df$nDNA != 1,]
to.fit$y = to.fit$meanmean/to.fit$nDNA
to.fit$x = 1- to.fit$mu*to.fit$nDNA**0.5
summary(lm(y ~ x, data=to.fit))
