library(ggplot2)
library(ggpubr)

for(scale in c("0.000", "0.500")) {
  for(penalty in c("0.000", "1.000")) {
    
    # read in all parallelised outputs for a given scale/penalty combination
    df = data.frame()
    for(expt in 0:6) {
      tdf = read.csv(paste0("inherit-mean-out-", expt, "-", scale, "-", penalty, ".csv"))
      df = rbind(df, tdf)
    }
    df = df[!is.na(df$expt),]
    
    my.trans = function(lambda) {
      return( ifelse(lambda == 0, log2(0.002), log2(lambda)))
    }
    
    ##### first look at the long-time behaviour (**NOT guaranteed to be steady state!)
    long.df = df[df$DUI == 0 & df$t == 1950,]
      
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
      scale_x_continuous(trans = "log", labels = scales::label_number(accuracy = 1e-4)) + 
      scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_tile() + facet_grid(mu ~ env)
    
    long.df = df[df$DUI == 1 & df$t == 1950,]
    
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
      scale_x_continuous(trans = "log", labels = scales::label_number(accuracy = 1e-4)) + 
      scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_tile() + facet_grid(mu ~ env)

    print(paste0("Outputting ", scale, " ", penalty))
    
    sf = 2
    png(paste0("panels-", scale, "-", penalty, ".png"), width=800*sf, height=1000*sf, res=72*sf)
    print( ggarrange(g.dui0, g.dui1, nrow=2, ncol=1, labels = c("A. DUI = 0", "B. DUI = 1")) )
    dev.off()
    
    #### now look at the time series to see if steady state looks reasonable
    t.df = df[df$DUI == 0 & df$t < 1950,]
    g.traces = ggplot(t.df, aes(x=t, y=mean.f/nDNA, color=factor(leakage+nDNA))) + geom_point(size=0.5) +
      facet_grid(mu ~ env) + theme(legend.position = "none")
    
    sf = 2
    png(paste0("panels-all-", scale, "-", penalty, ".png"), width=1000*sf, height=1000*sf, res=72*sf)
    print( ggarrange( ggarrange(g.dui0+theme(legend.position="none"), g.dui1+theme(legend.position="none"), nrow=1, ncol=2, labels = c("A. DUI = 0", "B. DUI = 1")),
               g.traces, nrow=2, ncol=1, labels=c("", "C. Time series") ) )
    dev.off()
    
  }
}

