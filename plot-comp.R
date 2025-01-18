library(ggplot2)
df = read.csv("inherit-comp-100-0.000-0.000-0.00000-1.000.csv")
ggplot(df, aes(x=t, y=mean.leak, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) 
ggplot(df, aes(x=t, y=mean.NDNA, color=factor(expt))) + geom_line() + facet_grid(mu ~ env) + xlim(0,5)
