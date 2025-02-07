library(ggplot2)
library(viridis)
library(dplyr)
library(ggpubr)
library(scales)

df.3 = read.csv("inherit-comp-1000-0.900-0.000-0.00000-1.000-1.000-0.csv")
df.4 = read.csv("inherit-comp-single-1000-0.900-0.000-0.00000-1.000-1.000-0.csv")

# previous plots for exploring simulation outcomes
if(FALSE) {
  ggplot(df, aes(x=t, y=mean.f, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) + theme(legend.position="none")
  ggplot(df, aes(x=t, y=mean.leak, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) 
  ggplot(df, aes(x=t, y=mean.NDNA, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) 
  
  ggplot(df[df$mu==1e-6 & df$env==16,], aes(x=t, y=mean.f, color=factor(expt))) + geom_line() + 
    facet_grid(mu ~ env) + xlim(0,100)
  
  ggplot(df, aes(x=t, y=mean.NDNA, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) + xlim(0,100)
  
  ggplot(df[df$t==499,], aes(x=log10(mean.NDNA), fill=factor(mu))) + geom_boxplot() + facet_wrap(~env) 
  ggplot(df[df$t==499,], aes(x=mean.leak, fill=factor(mu))) + geom_boxplot() + facet_wrap(~env) 
  ggplot(df[df$t==499,], aes(x=mean.leak, fill=factor(env))) + geom_boxplot() + facet_wrap(~mu) 
}

# summarise data from single-switch experiment
means.1 <- df.4[df.4$t==300,] %>%
  group_by(across(env:mu)) %>%  # Group by columns a through h
  summarize(mm.NDNA = mean(mean.NDNA, na.rm = TRUE), 
            mm.leak = mean(mean.leak, na.rm = TRUE), 
            mm.f = mean(mean.f, na.rm = TRUE),
            mm.h = mean(mean.h, na.rm = TRUE),
            .groups = 'drop')
means.1$mu = scientific(means.1$mu, digits=1)
means.1$mu = factor(means.1$mu, levels=unique(means.1$mu))

# summarise data from continuous-switching experiment
means.2 <- df.3[df.3$t==300,] %>%
  group_by(across(env:mu)) %>%  # Group by columns a through h
  summarize(mm.NDNA = mean(mean.NDNA, na.rm = TRUE), 
            mm.leak = mean(mean.leak, na.rm = TRUE), 
            mm.f = mean(mean.f, na.rm = TRUE),
            mm.h = mean(mean.h, na.rm = TRUE),
            .groups = 'drop')
means.2$mu = scientific(means.2$mu, digits=1)
means.2$mu = factor(means.1$mu, levels=unique(means.2$mu))

# fitness, mean N, mean lambda plots
plot.f.1= ggplot(means.1, aes(x=factor(mu), y=env, fill=mm.f)) + geom_tile() + 
  scale_y_continuous(transform="log2") + 
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.75, 0), limits=c(0.8,1))  +
  labs(x="Mutation rate", y="Env change period", fill="Popn mean\nfitness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot.n.1 = ggplot(means.1, aes(x=factor(mu), y=env, fill=mm.NDNA)) + geom_tile() + 
  scale_y_continuous(transform="log2") + 
  scale_fill_viridis(option="inferno") +
  labs(x="Mutation rate", y="Env change period", fill="Popn mean\nploidy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot.lambda.1 = ggplot(means.1, aes(x=factor(mu), y=env, fill=mm.leak)) + geom_tile() + 
  scale_y_continuous(transform="log2")+ 
  scale_fill_viridis(option="inferno") +
  labs(x="Mutation rate", y="Env change period", fill="Popn mean\nleakage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# fitness, mean N, mean lambda plots
plot.f.2= ggplot(means.2, aes(x=factor(mu), y=env, fill=mm.f)) + geom_tile() + 
  scale_y_continuous(transform="log2") + 
  scale_fill_gradientn(colors = c("black", "blue", "white", "red"), values = c(1, 1-1e-3, 0.75, 0), limits=c(0.8,1))  +
  labs(x="Mutation rate", y="Env change period", fill="Popn mean\nfitness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot.n.2 = ggplot(means.2, aes(x=factor(mu), y=env, fill=mm.NDNA)) + geom_tile() + 
  scale_y_continuous(transform="log2") + 
  scale_fill_viridis(option="inferno") +
  labs(x="Mutation rate", y="Env change period", fill="Popn mean\nploidy") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot.lambda.2 = ggplot(means.2, aes(x=factor(mu), y=env, fill=mm.leak)) + geom_tile() + 
  scale_y_continuous(transform="log2")+ 
  scale_fill_viridis(option="inferno") +
  labs(x="Mutation rate", y="Env change period", fill="Popn mean\nleakage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggarrange(plot.n.1, plot.lambda.1,
          plot.n.2, plot.lambda.2, nrow=2, ncol=2, labels=c("A", "B", "C", "D"))

sf = 2
png("optimal-behaviour.png", width=600*sf, height=400*sf, res=72*sf)
ggarrange(plot.n.1, plot.lambda.1,
          plot.n.2, plot.lambda.2, nrow=2, ncol=2, labels=c("A", "B", "C", "D"))

#ggarrange(plot.f, plot.n, plot.lambda, nrow=1)
dev.off()

# linear model for influences
summary(lm(mean.NDNA ~ env + mu, data=df[df$t==499,]))
summary(lm(mean.leak ~ env + mu, data=df[df$t==499,]))
