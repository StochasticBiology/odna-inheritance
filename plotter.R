library(ggplot2)
library(ggpubr)

npop = 50

###### phase space behaviour
# wrapper script currently only supports loop for npop = 100.
# for npop = 100, comment l8 and use l9-10; otherwise comment l9-10 and use l8 
scale = "0.500"; penalty = "0.000"; {{
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
  
  g.dui0 = ggplot(mean.df, aes(x=leakage,y=nDNA,fill=meanmean/nDNA)) + 
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0)) +
    geom_tile() + facet_grid(mu ~ env)
  
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
  
  g.dui1 = ggplot(mean.df, aes(x=leakage,y=nDNA,fill=meanmean/nDNA)) + 
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0)) +
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

##### "baseline" behaviour -- no mutation or DUI or env change

g.baseline = list()
for(scale in c("0.000", "0.500")) {
  for(penalty in c("0.000", "1.000")) {
    
    # read in all parallelised outputs for a given scale/penalty combination
    df = data.frame()
    for(expt in 0:1) {
      tdf = read.csv(paste0("inherit-baseline-mean-out-", npop, "-", expt, "-", scale, "-", penalty, ".csv"))
      df = rbind(df, tdf)
    }
    df = df[!is.na(df$expt),]
    
    long.df = df[df$DUI == 0 & df$t == 4950,]
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
    g.baseline[[length(g.baseline)+1]] = ggplot(mean.df, aes(x=mu, y=nDNA, fill=meanmean/nDNA)) + geom_tile() +
      scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) +
      scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0))
    
  }
}


ggarrange(plotlist = g.baseline, labels = c("A. 0, 0", "B. 0, 1", "C. 0.5, 0", "D. 0.5, 1"))

########### trying to collapse "baseline" behaviour into a trend

stat.df = mean.df
stat.df$response = mean.df$meanmean/mean.df$nDNA

my.lm = lm(log(response/(1-(response*0.9999))) ~ log(mu) + log(nDNA), 
           data=stat.df[stat.df$env==0,])
summary(my.lm)
#plot(my.lm)

ggplot(stat.df[stat.df$env==0,], aes(x=-1.11-1.18*log(mu)-1.81*log(nDNA), y=log(response/(1-response*0.9999)), color=factor(nDNA))) + geom_line() 

ggplot(stat.df[stat.df$env==0,], aes(x=-1.11-1.18*log(mu)-1.81*log(nDNA), y=log(response/(1-response)), color=factor(nDNA))) + geom_line() 

my.lm = lm(log(response/(1-(response*0.9999))) ~ log(mu) + log(nDNA), 
           data=stat.df[stat.df$env==0 & stat.df$nDNA <= 80,])
summary(my.lm)
#plot(my.lm)

ggplot(stat.df[stat.df$env==0,], aes(x=-1*log(mu)-2*log(nDNA), y=log(response/(1-response*0.9999)), color=factor(nDNA))) + geom_line() 
ggplot(stat.df[stat.df$env==0,], aes(x=-1*log(mu)-2*log(nDNA), y=log(response/(1-response*0.9999)), color=factor(nDNA))) + geom_line() 

baseline.fit.g = ggplot(stat.df[stat.df$env==0,], aes(x=-1.11-1.18*log(mu)-1.81*log(nDNA), y=log(response/(1-response)), color=factor(nDNA))) + geom_line() 

sf = 2
png(paste0("baseline-standard-", npop, ".png"), width=1000*sf, height=500*sf, res=72*sf)
print(ggarrange(g.baseline[[3]] + xlab("mu, 1e-6 = 0"),
                baseline.fit.g,
                labels = c("A. Baseline behaviour after 5000 gen", "B. Fit to logit-transformed heteroplasmy")))
dev.off()


########### zoomed-in region for simple case

# read in all parallelised outputs for a given scale/penalty combination
scale = "0.500"
penalty = "0.000"
df = data.frame()
for(expt in c(3,4,5,7,8,9,10)) {
  tdf = read.csv(paste0("oldzoom/inherit-zoom-mean-out-", expt, "-", scale, "-", penalty, ".csv"))
  df = rbind(df, tdf)
}
df = df[!is.na(df$expt),]
env.names = c("", "", "1", "10", "100", "", "5", "20", "50", "200")
df$period = as.numeric(env.names[df$env])

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

######## baseline behaviour for single oDNA per cell

scale = "0.500"; penalty = "0.000"
npop = 100; expt = 0

# read in all parallelised outputs for a given scale/penalty combination
df = data.frame()
tdf = read.csv(paste0("inherit-baseline-single-mean-out-", npop, "-", expt, "-", scale, "-", penalty, ".csv"))
df = rbind(df, tdf)
df = df[!is.na(df$expt),]

long.df = df[df$DUI == 0 & df$t == 4950,]
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
ggplot(mean.df, aes(x=mu, y=nDNA, fill=meanmean/nDNA)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0))

