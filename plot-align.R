library(reshape2)
library(ggplot2)
library(dplyr)

##### ICs matching BGP
# environmental period is longer here

# time series for zero mu
df = data.frame()
for(env in c(0, (1:10)*10, 150, 200)) {
  tdf = read.csv(paste0("inherit-belen-mean-out-100-", env, "-0.500-0.000.csv"))
  df = rbind(df, tdf)
}

ggplot(df[df$DUI==0 & df$mu==0,], aes(x=t, y=mean.f/nDNA, color=factor(nDNA))) + 
  geom_point(size=0.2, alpha=0.1) + facet_grid(env ~ leakage)

mdf = data.frame()
for(env in c(0, (1:10)*10, 150, 200)) {
  tdf = read.csv(paste0("inherit-belen-change-out-100-", env, "-0.500-0.000.csv"))
  mdf = rbind(mdf, tdf)
}

means.df = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))

ggplot(means.df, aes(x=log10(leakage), y=log10(nDNA), fill=mean_f)) + geom_tile() +
  facet_grid(env ~ mu)

#### better ICs
# environmental period is shorter here

# time series for zero mu
df = data.frame()
for(env in c(0, (1:10), 15, 20)) {
  tdf = read.csv(paste0("inherit-belen-mean-out-100-1-", env, "-0.500-0.000.csv"))
  df = rbind(df, tdf)
}

ggplot(df[df$DUI==0 & df$mu==0,], aes(x=t, y=mean.f/nDNA, color=factor(nDNA))) + 
  geom_point(size=0.2, alpha=0.1) + facet_grid(env ~ leakage)

# rough mean plots
mdf = data.frame()
for(env in c(0, (1:10), 15, 20)) {
  tdf = read.csv(paste0("inherit-belen-change-out-100-1-", env, "-0.500-0.000.csv"))
  mdf = rbind(mdf, tdf)
}

means.df = mdf %>% 
  group_by(leakage, DUI, nDNA, mu, env) %>%
  summarize(mean_f = mean(end.mean.f/nDNA))

ggarrange(
ggplot(means.df[means.df$DUI == 0,], aes(x=log10(leakage), y=log10(nDNA), fill=mean_f)) + geom_tile() +
  facet_grid(env ~ mu),
ggplot(means.df[means.df$DUI == 1,], aes(x=log10(leakage), y=log10(nDNA), fill=mean_f)) + geom_tile() +
  facet_grid(env ~ mu)
)

# empirical fit
sub = means.df[means.df$DUI == 0 & means.df$env == 0,]
ggplot(sub, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()

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


