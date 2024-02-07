library(ggplot2)
library(ggpubr)

# read in all parallelised outputs for a given scale/penalty combination
df = data.frame()
for(expt in 0:6) {
  tdf = read.csv(paste0("inherit-mean-out-", expt, "-0.500-0.000.csv"))
  df = rbind(df, tdf)
}
df = df[!is.na(df$expt),]

##### first look at the long-time behaviour (**NOT guaranteed to be steady state!)
long.df = df[df$DUI == 0 & df$t == 1950,]

mean.df = long.df[long.df$expt==0,]
mean.df$meanmean = 0
mean.df$meanvar = 0
for(i in 1:nrow(mean.df)) {
  mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]])
  mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                            long.df$nDNA == mean.df$nDNA[i] &
                                            long.df$mu == mean.df$mu[i] &
                                            long.df$leakage == mean.df$leakage[i]])
}

g.dui0 = ggplot(mean.df, aes(x=log10(leakage+0.005),y=log10(nDNA),fill=meanmean/nDNA)) + 
  geom_tile() + facet_grid(mu ~ env)

long.df = df[df$DUI == 1 & df$t == 1950,]

mean.df = long.df[long.df$expt==0,]
mean.df$meanmean = 0
mean.df$meanvar = 0
for(i in 1:nrow(mean.df)) {
  mean.df$meanmean[i] = mean(long.df$mean.f[long.df$env == mean.df$env[i] &
                                              long.df$nDNA == mean.df$nDNA[i] &
                                              long.df$mu == mean.df$mu[i] &
                                              long.df$leakage == mean.df$leakage[i]])
  mean.df$meanvar[i] = mean(long.df$var.f[long.df$env == mean.df$env[i] &
                                            long.df$nDNA == mean.df$nDNA[i] &
                                            long.df$mu == mean.df$mu[i] &
                                            long.df$leakage == mean.df$leakage[i]])
}

g.dui1 = ggplot(mean.df, aes(x=log10(leakage+0.005),y=log10(nDNA),fill=meanmean/nDNA)) + 
  geom_tile() + facet_grid(mu ~ env)

sf = 2
png("panels.png", width=800*sf, height=1000*sf, res=72*sf)
ggarrange(g.dui0, g.dui1, nrow=2, ncol=1)
dev.off()

#### now look at the time series to see if steady state looks reasonable
t.df = df[df$DUI == 0 & df$t < 1950,]
ggplot(t.df, aes(x=t, y=mean.f/nDNA, color=factor(leakage+nDNA))) + geom_point() +
  facet_grid(mu ~ env) + theme(legend.position = "none")
