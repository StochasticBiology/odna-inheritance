library(reshape2)
library(ggplot2)
library(dplyr)
library(metR)
library(ggpubr)

# explore effects of different model protocols 

mdf = data.frame()
for(env in c(0, 10)) {
  for(ics in c(0,1)) {
    for(leak in c(0,1)) {
      for(reamp in c(0,1)) {
        tdf = read.csv(paste0("inherit-belen-old-change-out-100-", ics, "-", env, "-0.500-0.000-", leak, "-", reamp, ".csv"))
        mdf = rbind(mdf, tdf)
      }
    }
  }
}

means.df = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env, ICs, det.leak, det.reamp) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))

expt.label = function(ICs, det.leak, det.reamp) {
  ics = c("hetero", "homo")[ICs]
  leak = c("s.leak", "d.leak")[det.leak]
  reamp = c("s.ramp", "d.ramp")[det.reamp]
  return(paste0(ics,leak,reamp, collapse=","))
}

means.df$label = 100*means.df$ICs+10*means.df$det.leak+means.df$det.reamp

sub = means.df[means.df$DUI==0 & means.df$env==0 & means.df$mu ==0, ]
g.env.0 = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)
sub = means.df[means.df$DUI==1 & means.df$env==0 & means.df$mu ==0, ]
g.env.0.dui = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)


sub = means.df[means.df$DUI==0 & means.df$env==10 & means.df$mu ==0, ]
g.env.10 = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)
sub = means.df[means.df$DUI==1 & means.df$env==10 & means.df$mu ==0, ]
g.env.10.dui = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)

ggarrange(g.env.0, g.env.10, 
          g.env.0.dui, g.env.10.dui,
          nrow=2, ncol=2, 
          labels = c("A. env=0, DUI=0", "B. env=10, DUI=0", "C. env=0, DUI=1", "D. env=10, DUI=1"))
# leading observations: 
# deterministic leakage model means high leakage always challenges adaptation (worse values with tens)
# heteroplasmic ICs maintain heterozygosity over longer time period (better values without hundreds for env = 10)

ggplot(means.df[means.df$env==10,]) +
  geom_point(aes(x=leakage*(1+0.1*DUI), y=mean_f, color=DUI, fill=nDNA)) +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  facet_wrap(~ label)
# DUI gives a modest advantage

#### finer steps through parameter space

## default experiment NB filename system updates
mdf = data.frame()
for(env in c(0, 2**(1:7))) {
  tdf = read.csv(paste0("inherit-template-change-out-100-1-", env, "-0.500-0.000-0-0-0.00000.csv"))
  mdf = rbind(mdf, tdf)
}

means.df.100 = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))

mdf = data.frame()
for(env in c(0, 2**(1:7))) {
  tdf = read.csv(paste0("inherit-template-change-out-100-1-", env, "-0.500-0.000-0-0-0.00100.csv"))
  mdf = rbind(mdf, tdf)
}

means.df.100.temp.1 = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))

mdf = data.frame()
for(env in c(0, 2**(1:7))) {
  tdf = read.csv(paste0("inherit-template-change-out-100-1-", env, "-0.500-0.000-0-0-0.10000.csv"))
  mdf = rbind(mdf, tdf)
}

means.df.100.temp.2 = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))

g.template = ggarrange(
  ggplot(means.df.100[means.df.100$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu),
  ggplot(means.df.100.temp.1[means.df.100.temp.1$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu),
  ggplot(means.df.100.temp.2[means.df.100.temp.2$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu),
  nrow=1
  
)

g.100 = ggarrange(
  ggplot(means.df.100[means.df.100$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu),
  ggplot(means.df.100[means.df.100$DUI == 1,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu)
)

## smaller population size

mdf = data.frame()
for(env in c(0, 2**(1:7))) {
  tdf = read.csv(paste0("inherit-belen-change-out-50-1-", env, "-0.500-0.000-0-0.csv"))
  mdf = rbind(mdf, tdf)
}

means.df.50 = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))


g.50 = ggarrange(
  ggplot(means.df.50[means.df.50$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu),
  ggplot(means.df.50[means.df.50$DUI == 1,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu)
)

## larger population size

mdf = data.frame()
for(env in c(0, 2**(1:7))) {
  tdf = read.csv(paste0("inherit-belen-change-out-200-1-", env, "-0.500-0.000-0-0.csv"))
  mdf = rbind(mdf, tdf)
}

means.df.200 = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))


g.200 = ggarrange(
  ggplot(means.df.200[means.df.200$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu),
  ggplot(means.df.200[means.df.200$DUI == 1,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu)
)

## heteroplasmy penalty

mdf = data.frame()
for(env in c(0, 2**(1:7))) {
  tdf = read.csv(paste0("inherit-belen-change-out-100-1-", env, "-0.500-0.250-0-0.csv"))
  mdf = rbind(mdf, tdf)
}

means.df.hetpen = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))


g.hetpen = ggarrange(
  ggplot(means.df.hetpen[means.df.hetpen$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu),
  ggplot(means.df.hetpen[means.df.hetpen$DUI == 1,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
    facet_grid(env ~ mu)
)

sf = 2
png("all-trellis.png", width=1400*sf, height=1200*sf, res=72*sf)
ggarrange(g.100, g.50, g.200, g.hetpen, ncol=2, nrow=2,
          labels =c("A. 100", "B. 50", "C. 200", "D. 100 hetpen"))
dev.off()

##### empirical fitting

means.df = means.df.100

model.R2 = function(x, y) {
  refs = which(!is.na(x) & !is.na(y))
  x = x[refs]
  y = y[refs]
  SSR = sum((x-y)**2)
  SST = sum((y-mean(y))**2)
  return(1-SSR/SST)
}

# empirical fit in constant environments
sub = means.df[means.df$DUI == 0 & means.df$env == 0 & means.df$leakage == 0,]
ggplot(sub, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()

sub = means.df[means.df$DUI == 0 & means.df$env == 0,]
ggplot(sub, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()
model.R2(sub$mean_f, 1-sub$mu*sqrt(sub$nDNA)/(1-sub$leakage)**2)

# remove the highest-leakage case, which is an outlier to this empirical fit
sub2 = sub[sub$leakage < 0.5,]
ggplot(sub2, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()
model.R2(sub2$mean_f, 1-sub2$mu*sqrt(sub2$nDNA)/(1-sub2$leakage)**2)

png("plot-align-empirical.png", width=400*sf, height=400*sf, res=72*sf)
ggplot(sub2, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()
dev.off()

### empirical where is best?
sf = 2
for(expt in c("normal", "hetpen", "template")) {
  if(expt == "normal") { means.df = means.df.100 }
  if(expt == "hetpen") { means.df = means.df.hetpen }
  if(expt == "template") { means.df = means.df.100.temp.2 }
  
  sub = means.df[means.df$DUI == 0 & !is.na(means.df$mean_f),]
  bests = data.frame()
  for(this.mu in unique(sub$mu)) {
    for(this.env in unique(sub$env)) {
      subsub = sub[sub$mu==this.mu & sub$env==this.env,]
      bestfitness = subsub[which(subsub$mean_f==max(subsub$mean_f)),]
      goodfitness = subsub[which(subsub$mean_f>=0.99*max(subsub$mean_f)),]
      bests = rbind(bests, data.frame(mu=this.mu, env=this.env, 
                                      bestn = mean(bestfitness$nDNA), bestlambda=mean(bestfitness$leakage),
                                      goodn = mean(goodfitness$nDNA), goodlambda=mean(goodfitness$leakage)))
    }
  }
  bests$env[bests$env==0] = 1
  bests$mu[bests$mu==0] = 1e-7
  g.n = ggplot(bests, aes(x=mu, y=env, fill=bestn)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  g.lambda = ggplot(bests, aes(x=mu, y=env, fill=bestlambda)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  g.n.good = ggplot(bests, aes(x=mu, y=env, fill=goodn)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  g.lambda.good = ggplot(bests, aes(x=mu, y=env, fill=goodlambda)) + geom_tile() +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  ggarrange(g.n, g.lambda,
            g.n.good, g.lambda.good)
  
  # LOESS fits to smooth contour plots
  mod.n = loess(goodn ~ mu + env, data = bests)
  to.plot.n = bests
  to.plot.n$goodn = predict(mod.n, newdata = to.plot.n, se = FALSE)
  mod.lambda = loess(goodlambda ~ mu + env, data = bests)
  to.plot.lambda = bests
  to.plot.lambda$goodlambda = predict(mod.lambda, newdata = to.plot.lambda, se = FALSE)
  
  # simulation data and smoother contours
  g.smooth.n = ggplot() +
    geom_tile(data=bests, aes(x=mu, y=env, fill=goodn)) +
    geom_contour(data=to.plot.n, aes(x=mu, y=env, z=goodn), 
                 breaks=c((1:10)*20), color="white", alpha=0.5) +
    geom_text_contour(data=to.plot.n, aes(x=mu, y=env, z=goodn), 
                      breaks=c((1:10)*20), color="white", skip=0, min.size = 0, alpha=0.5) +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  g.smooth.lambda = ggplot() +
    geom_tile(data=bests, aes(x=mu, y=env, fill=goodlambda)) +
    geom_contour(data=to.plot.lambda, aes(x=mu, y=env, z=goodlambda), 
                 breaks=c((0:5)*0.1), color="white", alpha=0.5) +
    geom_text_contour(data=to.plot.lambda, aes(x=mu, y=env, z=goodlambda), 
                      breaks=c((0:5)*0.1), color="white", skip=0, min.size = 0, alpha=0.5) +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  png(paste0("plot-align-optim-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(ggarrange(g.smooth.n, g.smooth.lambda))
  dev.off()
}

##### comparison of BGP and aligned dynamics -- for debugging only

# read a particular BGP example, produced by Inherit_comparison.py

dft = data.frame()
env = 8; mu = 0.01;
for(seed in c(0, 1, 2, 3, 4, 5, 6, 200, 500, 666)) {
  for(Ns in c(10*2**(0:5))) {
    for(lambda in c(0,"0.0005",0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.13,0.26,0.51)) {
      tdf = read.csv(paste0("f_env_", env, "_lambda_", lambda, "_mu_", mu, "_N_", Ns, "_seed_", seed, ".csv"))
      dft = rbind(dft, data.frame(meanf=mean(tdf[-1,499]),lambda=lambda,Ns=Ns,seed=seed))
    }
  }
}

dft$lambda = as.numeric(dft$lambda)
means.df = dft %>% 
  group_by(lambda, Ns) %>%
  summarize(mean_f = mean(meanf/Ns))

g.igj = ggplot(means.df.100[means.df.100$DUI == 0 & means.df.100$env==8 & means.df.100$mu == 0.01 & means.df.100$nDNA < 500,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) 

g.bgp = ggplot(means.df, aes(x=lambda, y=Ns, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) 

ggarrange(g.igj, g.bgp, labels=c("IGJ", "BGP"))

