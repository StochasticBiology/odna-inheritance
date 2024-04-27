library(reshape2)
library(ggplot2)
library(dplyr)

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


sub = means.df[means.df$DUI==0 & means.df$env==10 & means.df$mu ==0, ]
g.env.10 = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)

ggarrange(g.env.0, g.env.10, nrow=1)
# leading observations: 
# deterministic leakage model means high leakage always challenges adaptation (worse values with tens)
# heteroplasmic ICs maintain heterozygosity over longer time period (better values without hundreds for env = 10)


#### finer steps through parameter space

mdf = data.frame()
for(env in c(0, 2**(1:7))) {
  tdf = read.csv(paste0("inherit-belen-change-out-100-1-", env, "-0.500-0.000-0.csv"))
  mdf = rbind(mdf, tdf)
}

means.df = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))

### remaining question -- when is leakage bad in const envs??

ggarrange(
ggplot(means.df[means.df$DUI == 0,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
  facet_grid(env ~ mu),
ggplot(means.df[means.df$DUI == 1,], aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
  facet_grid(env ~ mu)
)

# empirical fit
sub = means.df[means.df$DUI == 0 & means.df$env == 0 & means.df$leakage == 0,]
ggplot(sub, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()

sub = means.df[means.df$DUI == 0 & means.df$env == 0,]
ggplot(sub, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(nDNA))) + 
  geom_point() + geom_abline()
cor(sub$mean_f, 1-sub$mu*sqrt(sub$nDNA)/(1-sub$leakage)**2, use="complete.obs" )**2



##### comparison of BGP and aligned dynamics -- for debugging

# read a particular BGP example

dft = data.frame()
env = 50; lambda = 0.01; Ns = 100; mu = 1e-6
for(seed in c(0, 1, 2, 3, 4, 5, 6, 200, 500, 666)) {
  tdf = read.csv(paste0("f_env_", env, "_lambda_", lambda, "_mu_", mu, "_N_", Ns, "_seed_", seed, ".csv"))
  tdft = as.data.frame(t(tdf)[-1,])
  colnames(tdft) = as.character(1:100)
  tdft$mean = rowMeans(tdft)
  tdft$t = 0:499
  
  dft = rbind(dft, tdft)
}

long_data <- melt(
  data = dft,
  id.vars = "t",
  variable.name = "variable",
  value.name = "value"
)

g.belen = ggplot(dft, aes(x=t, y=mean, color=factor(seed))) + geom_line()

# read an aligned example -- nb this particular param set may not appear in the set of simulated outputs

dfi = read.csv("inherit-belen-mean-out-100-50-0.500-0.000.csv")
dfi = dfi[dfi$DUI == 0, dfi$lambda = lambda, dfi$nDNA == Ns, dfi$mu == mu,]
g.igj = ggplot(dfi, aes(x=t, y=mean.f/nDNA, color=factor(nDNA))) + geom_line()

# compare
ggarrange(g.belen, g.igj, nrow=2)


