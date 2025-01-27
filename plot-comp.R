library(ggplot2)
df.1 = read.csv("inherit-comp-1000-0.000-0.000-0.00000-1.000-0.csv")
df.2 = read.csv("inherit-comp-1000-0.500-0.000-0.00000-1.000-0.csv")
df.3 = read.csv("inherit-comp-1000-0.900-0.000-0.00000-1.000-0.csv")

df = df.3

ggplot(df, aes(x=t, y=mean.f, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) + theme(legend.position="none")
ggplot(df, aes(x=t, y=mean.leak, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) 
ggplot(df, aes(x=t, y=mean.NDNA, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) 

ggplot(df[df$mu==1e-6 & df$env==16,], aes(x=t, y=mean.f, color=factor(expt))) + geom_line() + 
  facet_grid(mu ~ env) + xlim(0,100)

ggplot(df, aes(x=t, y=mean.NDNA, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) + xlim(0,100)

ggplot(df[df$t==499,], aes(x=log10(mean.NDNA), fill=factor(mu))) + geom_boxplot() + facet_wrap(~env) 
ggplot(df[df$t==499,], aes(x=mean.leak, fill=factor(mu))) + geom_boxplot() + facet_wrap(~env) 
ggplot(df[df$t==499,], aes(x=mean.leak, fill=factor(env))) + geom_boxplot() + facet_wrap(~mu) 

means.1 <- df[df$t==499,] %>%
  group_by(across(env:mu)) %>%  # Group by columns a through h
  summarize(mm.NDNA = mean(mean.NDNA, na.rm = TRUE), 
            mm.leak = mean(mean.leak, na.rm = TRUE), 
            mm.f = mean(mean.f, na.rm = TRUE),
            mm.h = mean(mean.h, na.rm = TRUE),
            .groups = 'drop')

ggplot(means.1, aes(x=log10(mu), y=env, fill=mm.NDNA)) + geom_tile() + 
  scale_y_continuous(transform="log2")
ggplot(means.1, aes(x=log10(mu), y=env, fill=mm.leak)) + geom_tile() + 
  scale_y_continuous(transform="log2")

summary(lm(mean.NDNA ~ env + mu, data=df[df$t==499,]))
summary(lm(mean.leak ~ env + mu, data=df[df$t==499,]))
