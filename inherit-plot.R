library(ggplot2)
library(grid)
library(gridExtra)

dfp = read.csv("inherit-out.csv")

exptstring = c("Neutral environment", "Environment favours A", "Environment favours B", "Rapid A-B switching env", "Slower A-B switching env")

for(env in c(0,1,3,4)) {
df = dfp[dfp$t %in% c(10, 100, 500, 900) & dfp$env == env,]
#df = dfp
df$propA = df$a/(df$a+df$b)
sub = df[df$DUI == 0 & df$LEAKAGE == 0,]
#g1a = ggplot(sub) + geom_line(aes(x=t,y=propA,color=factor(expt)))# + geom_line(aes(x=t,y=c/100,color=factor(expt)))
#g1 = ggplot(sub, aes(x=c,y=propA)) + geom_bin2d() + xlim(c(0,100)) + ylim(c(-0.5,1.5)) + facet_wrap(~t)
g1 = ggplot(sub, aes(x=c/100,y=propA,color=factor(MU))) + geom_jitter(width=0.01, height=0.01, alpha=0.2) + ylim(c(-0.1,1.1)) + facet_wrap(~t) + labs(color='0 DUI 0 leak') 
sub = df[df$DUI == 1 & df$LEAKAGE == 0,]
g2 = ggplot(sub, aes(x=c/100,y=propA,color=factor(MU))) + geom_jitter(width=0.01, height=0.01, alpha=0.2) + ylim(c(-0.1,1.1)) + facet_wrap(~t) + labs(color='1 DUI 0 leak') 
#g2 = ggplot(sub, aes(x=c,y=propA)) + geom_bin2d() + xlim(c(0,100)) + ylim(c(-0.5,1.5)) + facet_wrap(~t)
sub = df[df$DUI == 0 & df$LEAKAGE == 0.05,]
g3 = ggplot(sub, aes(x=c/100,y=propA,color=factor(MU))) + geom_jitter(width=0.01, height=0.01, alpha=0.2)  + ylim(c(-0.1,1.1)) + facet_wrap(~t) + labs(color='0 DUI 1 leak') 
#g3 = ggplot(sub, aes(x=c,y=a)) + geom_bin2d() + xlim(c(0,100)) + ylim(c(0,100)) + facet_wrap(~t)
sub = df[df$DUI == 1 & df$LEAKAGE == 0.05,]
g4 = ggplot(sub, aes(x=c/100,y=propA,color=factor(MU))) + geom_jitter(width=0.01, height=0.01, alpha=0.2)  +ylim(c(-0.1,1.1)) + facet_wrap(~t) + labs(color='1 DUI 1 leak') 
#g4 = ggplot(sub, aes(x=c,y=a)) + geom_bin2d() + xlim(c(0,100)) + ylim(c(0,100)) + facet_wrap(~t)

png(paste(c("inherit-env-", env, ".png"), collapse=""), width=800, height=600)
grid.arrange(g1,g2,g3,g4,nrow=2,top = textGrob(exptstring[env+1],gp=gpar(fontsize=20,font=3)))
dev.off()
}

for(env in c(0,1,3,4)) {
  df = dfp[dfp$env == env,]
 
  sub = df[df$DUI == 0 & df$LEAKAGE == 0,]
  g1 = ggplot(sub, aes(x=t,y=f,color=factor(expt),fill=factor(expt))) + geom_smooth(alpha=0.05) + theme(legend.position = "none") # labs(color='0 DUI 0 leak') 
  sub = df[df$DUI == 1 & df$LEAKAGE == 0,]
  g2 = ggplot(sub, aes(x=t,y=f,color=factor(expt),fill=factor(expt))) + geom_smooth(alpha=0.05) +  theme(legend.position = "none") # labs(color='1 DUI 0 leak') 
  #g2 = ggplot(sub, aes(x=c,y=propA)) + geom_bin2d() + xlim(c(0,100)) + ylim(c(-0.5,1.5)) + facet_wrap(~t)
  sub = df[df$DUI == 0 & df$LEAKAGE == 0.05,]
  g3 = ggplot(sub, aes(x=t,y=f,color=factor(expt),fill=factor(expt))) + geom_smooth(alpha=0.05) +  theme(legend.position = "none") # labs(color='0 DUI 1 leak') 
   #g3 = ggplot(sub, aes(x=c,y=a)) + geom_bin2d() + xlim(c(0,100)) + ylim(c(0,100)) + facet_wrap(~t)
  sub = df[df$DUI == 1 & df$LEAKAGE == 0.05,]
  g4 = ggplot(sub, aes(x=t,y=f,color=factor(expt),fill=factor(expt))) + geom_smooth(alpha=0.05) +  theme(legend.position = "none") # labs(color='1 DUI 1 leak') 
  #g4 = ggplot(sub, aes(x=c,y=a)) + geom_bin2d() + xlim(c(0,100)) + ylim(c(0,100)) + facet_wrap(~t)
  
  png(paste(c("inherit-trace-", env, ".png"), collapse=""), width=800, height=600)
  grid.arrange(g1,g2,g3,g4,nrow=2,top = textGrob(exptstring[env+1],gp=gpar(fontsize=20,font=3)))
  dev.off()
}

#grid.arrange(g1,g2,nrow=1)

