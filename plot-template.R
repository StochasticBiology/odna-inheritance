library(reshape2)
library(ggplot2)
library(viridis)
library(dplyr)
library(metR)
library(ggpubr)

# time series demonstrations

df1 = read.csv("inherit-template-full-out-10-1-0-0.500-0.000-0-0-0.00000.csv")
expt.1 = df1[df1$t < 50,]

df2 = read.csv("inherit-template-full-out-10-1-8-0.500-0.000-0-0-0.00000.csv")
expt.2 = df2[df2$t < 50,]

# function producting a plot summarising dynamics for a given subset of observations
dyn.sum.plot = function(plot.set, tau, tstr) {
   means.1 <- plot.set %>%
    group_by(across(Npop:t)) %>%  # Group by columns a through h
    summarize(mean_a = mean(a, na.rm = TRUE), 
              mean_b = mean(b, na.rm = TRUE), 
              mean_c = mean(c, na.rm = TRUE),
              mean_f = mean(f, na.rm = TRUE),
              .groups = 'drop')
  
 return( ggplot() +
           geom_rect(data=data.frame(xmin=0,ymin=0,xmax=tau,ymax=50), 
                     aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill="#8888FF", alpha=0.4) +
           geom_rect(data=data.frame(xmin=tau,ymin=0,xmax=35,ymax=50), 
                     aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill="#FF8888", alpha=0.4) +
    geom_line(data=means.1, aes(x=t, y=mean_a), color="blue", linewidth=2) +
    geom_line(data=means.1, aes(x=t, y=mean_b), color="red", linewidth=2) +
    geom_line(data=means.1, aes(x=t, y=mean_c), color="grey", linewidth=2) +
     
    geom_point(data=plot.set, aes(x=t, y=a), color="blue", size=2, alpha=0.5) +
    geom_point(data=plot.set, aes(x=t+0.3, y=b), color="red", size=2, alpha=0.5) +
    geom_point(data=plot.set, aes(x=t+0.6, y=c), color="grey", size=2, alpha=0.5) +
      geom_line(data=means.1, aes(x=t, y=mean_f), color="black", linewidth=2) +
      labs(x="t", y="Counts and fitness") +
      xlim(0,35) +
      ggtitle(tstr) +
    theme_light()
 )
  
}


# output a given set of example trajectories
sf =2
png("ex-time-series.png", width=600*sf, height=400*sf, res=72*sf)
ggarrange(
  dyn.sum.plot(expt.1[expt.1$mu==0 & expt.1$leakage==0.0,], 0, "μ=0, λ=0, τ=1"),
  dyn.sum.plot(expt.1[expt.1$mu==2e-4 & expt.1$leakage==0.25,], 0, "μ=2e-4, λ=0.25, τ=1"),
  dyn.sum.plot(expt.2[expt.2$mu==0 & expt.1$leakage==0.,], 8, "μ=0, λ=0, τ=8"),
  dyn.sum.plot(expt.2[expt.2$mu==2e-4 & expt.1$leakage==0.25,], 8, "μ=2e-4, λ=0.25, τ=8"),
  labels = c("B", "C", "D", "E" )
)
dev.off()


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

label.fn = function(dfrow) {
  if(dfrow["ICs"] == 1) { s1 = "pIC"} else { s1 = "hIC"}
  if(dfrow["det.leak"] == 1) { s2 = "dL" } else { s2 = "rL"}
  if(dfrow["det.reamp"] == 1) { s3 = "dR" } else { s3 = "rR" }
  if(dfrow["ICs"] == 1 & dfrow["det.leak"] == 0 & dfrow["det.reamp"] == 0) { s4 = "*"} else {s4 = ""}
  return(paste0(c(s1,s2,s3,s4),collapse=" "))
}

means.df$label = 100*means.df$ICs+10*means.df$det.leak+means.df$det.reamp
means.df$label = apply(means.df, 1, function(dfrow) label.fn(dfrow))

sub = means.df[means.df$DUI==0 & means.df$env==0 & means.df$mu ==0, ]
g.env.0 = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
  labs(x = "Leakage", y="N", fill = "Mean\nfitness") +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)
sub = means.df[means.df$DUI==1 & means.df$env==0 & means.df$mu ==0, ]
g.env.0.dui = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Leakage", y="N", fill = "Mean\nfitness") +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)


sub = means.df[means.df$DUI==0 & means.df$env==10 & means.df$mu ==0, ]
g.env.10 = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Leakage", y="N", fill = "Mean\nfitness") +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)
sub = means.df[means.df$DUI==1 & means.df$env==10 & means.df$mu ==0, ]
g.env.10.dui = ggplot(sub, aes(x=leakage, y=nDNA, fill=mean_f)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Leakage", y="N", fill = "Mean\nfitness") +
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)) +
  #  facet_wrap(~expt.label(ICs, det.leak, det.reamp))
  facet_wrap(~label)

sf = 2
png("si-different-models.png", width=600*sf, height=600*sf, res=72*sf)
ggarrange(g.env.0, g.env.10, 
          g.env.0.dui, g.env.10.dui,
          nrow=2, ncol=2, 
          #labels = c("A. env=0, DUI=0", "B. env=10, DUI=0", "C. env=0, DUI=1", "D. env=10, DUI=1"))
          labels = c("A", "B", "C", "D"))
dev.off()

# leading observations: 
# deterministic leakage model means high leakage always challenges adaptation (worse values with tens)
# heteroplasmic ICs maintain heterozygosity over longer time period (better values without hundreds for env = 10)

ggplot(means.df[means.df$env==10,]) +
  geom_point(aes(x=leakage*(1+0.1*DUI), y=mean_f, color=DUI, fill=nDNA)) +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  facet_wrap(~ label)
# DUI gives a modest advantage

############## main block for visualisations

#### finer steps through parameter space

# general function to read in data output
read.datafile = function(NPOP = 100, ICs = 1, scale = 0.5, 
                         penalty = 0, TEMPLATE = 0, FITNESSB = 1, FITNESSC = 1) {
  mdf = data.frame()
  for(env in c(0, 2**(1:7))) {
    tdf = read.csv(sprintf("inherit-template-change-out-%i-%i-%i-%.3f-%.3f-0-0-%.5f-%.2f-%.2f.csv", 
                           NPOP, ICs, env, scale, penalty, TEMPLATE, FITNESSB, FITNESSC))
    mdf = rbind(mdf, tdf)
  }
  return( mdf %>% 
    group_by(leakage, DUI, nDNA, mu, env) %>%
    summarize(mean_f = mean(end.mean.f/nDNA))
  )
}

## default experiment NB filename system updates

means.df.100 = read.datafile()
## templated repair rate 1
means.df.100.temp.1 = read.datafile(TEMPLATE = 0.001)
## templated repair rate 2
means.df.100.temp.2 = read.datafile(TEMPLATE = 0.1)

## smaller population size
means.df.50 = read.datafile(NPOP = 50)
## larger population size
means.df.200 = read.datafile(NPOP = 200)

## heteroplasmy penalty
means.df.hetpen = read.datafile(penalty = 0.25)

###############
###### this section added after peer review: looking at different within-cell fitnesses

## postivie/negative cellular fitness for B/C
means.df.fit.neg.b = read.datafile(FITNESSB = 0.91)
means.df.fit.pos.b = read.datafile(FITNESSB = 1.1)
means.df.fit.neg.c = read.datafile(FITNESSC = 0.91)
means.df.fit.pos.c = read.datafile(FITNESSC = 1.1)
means.df.fit.pos.c.template = read.datafile(TEMPLATE = 0.001, FITNESSC = 1.1)

mu.set = c(0,0.0001,0.01)
env.set = c(2, 8, 16, 64)
strip.pos = "top"
arrangement = facet_grid(env ~ mu)

g.format = list(geom_tile(aes(x=leakage, y=nDNA, fill=mean_f)),
                scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)),
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)),
  theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  labs(x = "Leakage", y = "N", fill = "Mean\nfitness"),
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.5, 0), limits=c(0,1)),
  arrangement)

pull.set = function(df, DUI = 0) {
  return(df[df$mu %in% mu.set & df$env %in% env.set & df$DUI == DUI,])
}

g.sub.fit.neg.b = ggplot(pull.set(means.df.fit.neg.b)) + g.format + ggtitle("Cell disadv to B")
g.sub.fit.pos.b = ggplot(pull.set(means.df.fit.pos.b)) + g.format + ggtitle("Cell adv to B")
g.sub.fit.neg.c = ggplot(pull.set(means.df.fit.neg.c)) + g.format + ggtitle("Cell disadv to C")
g.sub.fit.pos.c = ggplot(pull.set(means.df.fit.pos.c)) + g.format + ggtitle("Cell adv to C")

ggarrange(g.sub.fit.neg.b, g.sub.fit.pos.b,
          g.sub.fit.neg.c, g.sub.fit.pos.c)

###### end of new section
###############

g.all.100.0 = ggplot(means.df.100[means.df.100$DUI == 0,]) +g.format

wdf = means.df.100[means.df.100$DUI == 0,]
wdf$pred = -log(1+wdf$env) / (log(1+wdf$leakage*wdf$nDNA)) - log(1+wdf$mu)*log(1+wdf$nDNA)

ggplot(wdf, aes(x=leakage, y=nDNA, fill=pred)) + geom_tile() +
  scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-3, 1e-2, 1e-1, 1)) +
  scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(10, 20, 50, 100, 200, 500)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(env ~ mu)

for(expt.detail in c(0, 1)) {
  if(expt.detail == 0) {
    mu.set = c(0.01)
    env.set = c(64)
    strip.pos = "none"
    arrangement = NULL
  } else {
    mu.set = c(0,0.0001,0.01)
    env.set = c(2, 8, 16, 64)
    strip.pos = "top"
    arrangement = facet_grid(env ~ mu)
  }
  
  g.sub.100.0 = ggplot(pull.set(means.df.100)) + g.format + ggtitle("Standard")
  g.sub.100.1 = ggplot(pull.set(means.df.100, DUI=1)) + g.format + ggtitle("DUI")
  g.sub.100.temp.1.0 = ggplot(pull.set(means.df.100.temp.1)) + g.format + ggtitle("Templated repair")
  g.sub.100.temp.2.0 = ggplot(pull.set(means.df.100.temp.2)) + g.format + ggtitle("Templated repair alt")

  g.sub.100.temp.1.1 = ggplot(pull.set(means.df.100.temp.1, DUI=1)) + g.format + ggtitle("DUI, templated repair")

  g.sub.50.0 = ggplot(pull.set(means.df.50)) + g.format + ggtitle("Npop = 50")
  g.sub.50.1 = ggplot(pull.set(means.df.50, DUI=1)) + g.format + ggtitle("DUI, Npop = 50")
  g.sub.200.0 = ggplot(pull.set(means.df.200)) + g.format + ggtitle("Npop = 200")
  g.sub.200.1 = ggplot(pull.set(means.df.200, DUI=1)) + g.format + ggtitle("DUI, Npop = 200")
  
  g.sub.hetpen.0 = ggplot(pull.set(means.df.hetpen)) + g.format + ggtitle("h penalty")
  g.sub.hetpen.1 = ggplot(pull.set(means.df.hetpen, DUI=1)) + g.format + ggtitle("DUI, h penalty")
  
  g.sub.fit.neg.b = ggplot(pull.set(means.df.fit.neg.b)) + g.format + ggtitle("Cell disadv to B")
  g.sub.fit.pos.b = ggplot(pull.set(means.df.fit.pos.b)) + g.format + ggtitle("Cell adv to B")
  g.sub.fit.neg.c = ggplot(pull.set(means.df.fit.neg.c)) + g.format + ggtitle("Cell disadv to M")
  g.sub.fit.pos.c = ggplot(pull.set(means.df.fit.pos.c)) + g.format + ggtitle("Cell adv to M")
  g.sub.fit.pos.c.template = ggplot(pull.set(means.df.fit.pos.c.template)) + g.format + ggtitle("Cell adv to M, templated repair")
  
  nl = theme(legend.position = "none")
  sf = 2
  if(expt.detail == 0) {
    png("all-trellises-update.png", width=900*sf, height=900*sf, res=72*sf)
    print(ggarrange(g.sub.100.0+nl, g.sub.50.0+nl, g.sub.200.0+nl, 
                    g.sub.100.1+nl, g.sub.100.temp.1.0+nl, g.sub.hetpen.0+nl, 
                    g.sub.fit.neg.b+nl, g.sub.fit.pos.c+nl, g.sub.fit.pos.c.template+nl,
                    #labels=c("A. Default", "B. Small popn", "C. Large popn", "D. DUI", "E. Repair", "F. h penalty"),
                    labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                    ncol=3,nrow=3))
    dev.off()
  } else {
    png("all-trellises-update-full.png", width=1200*sf, height=1200*sf, res=72*sf)
    print(ggarrange(g.sub.100.0+nl, g.sub.50.0+nl, g.sub.200.0+nl, 
                    g.sub.100.1+nl, g.sub.100.temp.1.0+nl, g.sub.hetpen.0+nl, 
                    g.sub.hetpen.1+nl, 
                    g.sub.fit.neg.b+nl, g.sub.fit.pos.b+nl, g.sub.fit.neg.c+nl, g.sub.fit.pos.c+nl,
                    g.sub.fit.pos.c.template+nl,
#                    labels=c("A. Default", "B. Small popn", "C. Large popn", "D. DUI", "E. Repair", "F. h penalty"),
                    labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
                    ncol=4,nrow=3))
    dev.off()
  }
 

}

png("main-text-trellis.png", width=600*sf, height=600*sf, res=72*sf)
print(g.all.100.0 + labs(x="Leakage", y="N", fill="Mean\nfitness"))
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

roze.model = function(mu, leakage, nDNA, norg, sc, sm, p) {
  alpha = leakage
  B = nDNA
  N = 2*B
  n = norg
  r = 1 - (2*alpha*(1-alpha)*(1-1/N))/(1 - (1-1/B)*(1-1/N)*(1 - 2*alpha*(1-alpha)))
  Ne = n / (1 + (1-r)*(1 - 2*(1-1/N)*(1-1/B)))
  S = sc + (1-1/B)*(1-r)*(sm - sc*(1-1/N))
  u = (1-exp(-2*Ne*S*p))/(1-exp(-2*Ne*S))
  return(u)
}

# empirical fit in constant environments
sub = means.df[means.df$DUI == 0 & means.df$env == 0 & means.df$leakage == 0,]
sub$roze = roze.model(sub$mu, sub$leakage, sub$nDNA, 100, 1, 0.5, 0.5)
ggplot(sub, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()

sub = means.df[means.df$DUI == 0 & means.df$env == 0,]
ggplot(sub, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()
model.R2(sub$mean_f, 1-sub$mu*sqrt(sub$nDNA)/(1-sub$leakage)**2)

# remove the highest-leakage case, which is an outlier to this empirical fit
sub2 = sub[sub$leakage < 0.5,]
sub2$roze = roze.model(sub2$mu, sub2$leakage, sub2$nDNA, 100, 1, 0.5, 0.5)
ggplot(sub2, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()
model.R2(sub2$mean_f, 1-sub2$mu*sqrt(sub2$nDNA)/(1-sub2$leakage)**2)
# align with roze
model.R2(sub2$mean_f, 1-sub2$mu*sqrt(sub2$nDNA)/(sub2$lambda*(1-sub2$leakage))/(1-(sub2$lambda*(1-sub2$leakage))))

ggplot(sub2, aes(x=roze, y=mean_f, color=factor(leakage))) + 
  geom_point() + geom_abline()

png("plot-align-empirical.png", width=400*sf, height=300*sf, res=72*sf)
ggplot(sub2, aes(x=1-mu*sqrt(nDNA)/(1-leakage)**2, y=mean_f, color=factor(leakage))) + 
  theme_light() + labs(color="Leakage", x="Empirical prediction", y="Mean fitness in simulation") +
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
  mod.n = loess(goodn ~ mu + env, data = bests, span = 0.5)
  to.plot.n = bests
  to.plot.n$goodn = predict(mod.n, newdata = to.plot.n, se = FALSE)
  mod.lambda = loess(goodlambda ~ mu + env, data = bests)
  to.plot.lambda = bests
  to.plot.lambda$goodlambda = predict(mod.lambda, newdata = to.plot.lambda, se = FALSE)
  
  # simulation data and smoother contours
  g.smooth.n = ggplot() +
    geom_tile(data=bests, aes(x=mu, y=env, fill=goodn)) +
    geom_contour(data=to.plot.n, aes(x=mu, y=env, z=goodn), 
                 breaks=c((1:10)*20), color="black", alpha=0.5) +
    geom_text_contour(data=to.plot.n, aes(x=mu, y=env, z=goodn), 
                      breaks=c((1:10)*20), color="black", skip=0, min.size = 0, alpha=1) +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    scale_fill_continuous(limits=c(0,200)) +
    scale_fill_viridis(option="inferno") +
    theme_light() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Mutation rate", y = "Environmental\nchange period", fill="Mean\noptimal\nN")
  
  g.smooth.lambda = ggplot() +
    geom_tile(data=bests, aes(x=mu, y=env, fill=goodlambda)) +
    geom_contour(data=to.plot.lambda, aes(x=mu, y=env, z=goodlambda), 
                 breaks=c((0:5)*0.1), color="black", alpha=0.5) +
    geom_text_contour(data=to.plot.lambda, aes(x=mu, y=env, z=goodlambda), 
                      breaks=c((0:5)*0.1), color="black", skip=0, min.size = 0, alpha=1) +
    scale_x_continuous(trans = "log", labels = function(x) format(x, scientific = TRUE), breaks = c(1e-6, 1e-4, 1e-2, 1)) +
    scale_y_continuous(trans = "log", labels = scales::label_number(accuracy = 1), breaks=c(2**(1:7))) +
    scale_fill_continuous(limits=c(0,0.55)) +
    scale_fill_viridis(option="inferno") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(x = "Mutation rate", y = "Environmental\nchange period", fill="Mean\noptimal\nleakage")
  
  png(paste0("plot-align-optim-", expt, ".png", collapse=""), width=600*sf, height=300*sf, res=72*sf)
  print(ggarrange(g.smooth.n, g.smooth.lambda, labels=c("A", "B")))
  dev.off()
}

##### comparison of BGP and aligned dynamics -- for debugging only, not used in final article
if(FALSE) {
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
}
