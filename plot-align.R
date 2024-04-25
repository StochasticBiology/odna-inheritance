library(reshape2)
library(ggplot2)

dft = data.frame()
for(seed in c(0, 1, 2, 3, 4, 5, 6, 200, 500, 666)) {
  for(env in c(0,50,70,100,150)) {
    for(lambda in c(0, 0.01, 0.05, 0.25)){
      for(Ns in c(10, 50, 100, 200)) {
        for(mu in c(1e-6, 1e-3, 1e-2)) {
   tdf = read.csv(paste0("f_env_", env, "_lambda_", lambda, "_mu_", mu, "_N_", Ns, "_seed_", seed, ".csv"))
tdft = as.data.frame(t(tdf)[-1,])
colnames(tdft) = as.character(1:100)
tdft$mean = rowMeans(tdft)
tdft$t = 0:499
tdft$seed = seed
tdft$env = env
tdft$lambda = lambda
tdft$mu = mu
tdft$N = Ns
dft = rbind(dft, tdft)
        }
      }
    }
  }
}

#df = read.csv("f_env_50_lambda_0_mu_1e-06_N_10_seed_666.csv")
#nrow(df)

long_data <- melt(
  data = dft,
  id.vars = "t",
  variable.name = "variable",
  value.name = "value"
)

#ggplot(long_data, aes(x=t, y=value, color=factor(variable))) + geom_point() + theme(legend.position="none")

g.belen = ggplot(dft, aes(x=t, y=mean, color=factor(seed))) + geom_line()

##### matching BGP

df = data.frame()
for(env in c(0, (1:10)*10, 150, 200)) {
  tdf = read.csv(paste0("inherit-belen-mean-out-100-", env, "-0.500-0.000.csv"))
df = rbind(df, tdf)
}

ggplot(df[df$DUI==0 & df$mu==0,], aes(x=t, y=mean.f/nDNA, color=factor(nDNA))) + 
  geom_point(size=0.2, alpha=0.1) + facet_grid(env ~ leakage)

#### better ICs

df = data.frame()
for(env in c(0, (1:10), 15, 20)) {
  tdf = read.csv(paste0("inherit-belen-mean-out-100-1-", env, "-0.500-0.000.csv"))
  df = rbind(df, tdf)
}

ggplot(df[df$DUI==0 & df$mu==0,], aes(x=t, y=mean.f/nDNA, color=factor(nDNA))) + 
  geom_point(size=0.2, alpha=0.1) + facet_grid(env ~ leakage)

#####

dfi = read.csv("inherit-belen-mean-out-100-50-0.500-0.000.csv")
#dfi = read.csv("inherit-belen-mean-out-10-20-0.000-0.000.csv")
g.igj = ggplot(dfi, aes(x=t, y=mean.f/nDNA, color=factor(nDNA))) + geom_line()

ggarrange(g.belen, g.igj, nrow=2)

####

dffull = read.csv("inherit-belen-full-out-100-50-0.500-0.000.csv")
g.a.igj = ggplot(dffull, aes(x=t, y=a, color=factor(expt))) + geom_point()

tadf = read.csv(paste0("a_env_50_lambda_0.2_mu_0_N_100_seed_", 666, ".csv"))
tdaft = as.data.frame(t(tadf)[-1,])
colnames(tdaft) = as.character(1:100)
tdaft$t = 0:499
long_data_a <- melt(
  data = tdaft,
  id.vars = "t",
  variable.name = "variable",
  value.name = "value"
)
g.a.belen = ggplot(long_data_a, aes(x=t, y=value)) + geom_point()

ggarrange(g.a.igj, g.a.belen)

#######

dfss = read.csv("a_env_50_lambda_0.2_mu_0_N_100_seed_6.csv")
