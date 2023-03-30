data <- (read.table("/Users/batistapj/Desktop/miRNA_ecdf", header = T))
head(data)

plot(ecdf(data$l2fc_200), verticals=TRUE, do.p=FALSE, main="ECDF plot for both samples", xlab="Scores", ylab="Cumulative Percent",lty="dashed")
lines(ecdf(data$l2fc_nT_200), verticals=TRUE, do.p=FALSE, col.h="red", col.v="red",lty="dotted")
lines(ecdf(data$l2fc_3065), verticals=TRUE, do.p=FALSE, col.h="green", col.v="green",lty="dotted")
lines(ecdf(data$l2fc_142), verticals=TRUE, do.p=FALSE, col.h="orange", col.v="orange",lty="dotted")
lines(ecdf(data$l2fc_nonT_all), verticals=TRUE, do.p=FALSE, col.h="blue", col.v="blue",lty="dotted")
lines(ecdf(data$l2fc_all), verticals=TRUE, do.p=FALSE, col.h="red", col.v="red",lty="dotted")

boxplot(data$l2fc_200, data$l2fc_3065, data$l2fc_142, data$l2fc_nonT_all)

vio<- (read.table("/Users/batistapj/Desktop/violoin_miR.txt", header = T))
head(vio)
library(scales)
library(ggplot2)
v <- vio
head(v)
p <- ggplot(v, aes(type,l2fc))
p + geom_boxplot()
p + geom_violin() + scale_y_continuous(trans=log2_trans())
p + geom_violin() + coord_flip()
p + geom_violin(scale = "count")
p + geom_violin(scale = "width")
p + geom_violin(trim = FALSE)
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun=mean, geom="point", shape=23, size=2)
p + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun=median, geom="point", size=2, color="red")
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun.data=mean_sdl, geom="crossbar", width=0.1)
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")

p + geom_violin(aes(fill = factor(type)))

p + geom_boxplot(aes(fill = factor(type))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")

###final version

v<- (read.table("/Users/batistapj/Desktop/violoin_final_miR.txt", header = T))
head(v)

p <- ggplot(v, aes(type,l2fc))
p + geom_boxplot(aes(fill = factor(type)))
p + geom_violin(aes(fill = factor(type)))

###STAT test

test_Ctrl <- (data$l2fc_all)
test_miR200 <- (data$l2fc_200)
test_noT <- (data$l2fc_nonT_all)

plot(ecdf(test_Ctrl), xlim = range(c(x, x2)))
plot(ecdf(test_miR200), add = TRUE, lty = "dotted")
ks.test(test_Ctrl, test_miR200)
ks.test(test_miR200, test_Ctrl, alternative = "g")
ks.test(test_miR200, test_Ctrl, alternative = "l")
ks.test(test_miR200, test_Ctrl, alternative = "two.sided")

ks.test(test_noT, test_Ctrl, alternative = "g")

wilcox.test (test_miR200, test_Ctrl)

# test if x is stochastically larger than x2
x <- rnorm(50)
x2 <- rnorm(100, -1)
plot(ecdf(x), xlim = range(c(x, x2)))
plot(ecdf(x2), add = TRUE, lty = "dashed")
t.test(x, x2, alternative = "g")
wilcox.test(x, x2, alternative = "g")
ks.test(x, x2, alternative = "l")


