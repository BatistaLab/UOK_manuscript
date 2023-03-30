

Target <- (read.table("/Users/batistapj/Desktop/ecdf/target7.txt", header = T))
Nontarget<- (read.table("/Users/batistapj/Desktop/ecdf/nontarget7.txt", header = T))
tmore2<-(read.table("/Users/batistapj/Desktop/ecdf/tmore2.txt", header = T))

head(Target)
dim(Target)
head(Nontarget)
dim(Nontarget)

plot(ecdf(Target$log2FoldChange), verticals=TRUE, do.p=FALSE, main="ECDF plot for both samples", xlab="Scores", ylab="Cumulative Percent",lty="dashed")
lines(ecdf(Nontarget$log2FoldChange), verticals=TRUE, do.p=FALSE, col.h="red", col.v="red",lty="dotted")
lines(ecdf(tmore2$log2FoldChange), verticals=TRUE, do.p=FALSE, col.h="red", col.v="red",lty="dotted")

boxplot(Target$log2FoldChange, Nontarget$log2FoldChange, tmore2$log2FoldChange)

vio<- (read.table("/Users/batistapj/Documents/Collaborations/YTHDC2_project/Review_reanalysis/ecdf/violin.txt", header = T))
library(scales)

v <- vio
head(v)
p <- ggplot(v, aes(type,L2FC))
p + geom_boxplot()
p + geom_violin() + scale_y_continuous(trans=log2_trans())
p + geom_violin() + coord_flip()
p + geom_violin(scale = "count")
p + geom_violin(scale = "width")
p + geom_violin(trim = FALSE)
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
p + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun.y=median, geom="point", size=2, color="red")
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun.data=mean_sdl, geom="crossbar", width=0.1)
p + geom_violin(aes(fill = factor(type))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")

p + geom_boxplot(aes(fill = factor(type))) + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")


#have to figure out how to switch the analysis to e in reference to wild type

test1 <- (AGOrpm1$log2FoldChange)
test2 <- (AGOrpm0$V3)
wilcox.test (test1, test2)
