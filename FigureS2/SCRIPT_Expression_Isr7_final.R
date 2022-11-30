### Marion Müller

library(edgeR)
library(dplyr)
library(ggplot2)

count <- read.table("Quant_Isr7_3rep.txt", header=TRUE, sep="\t", row.names=1)
head(count)

count_sub <- count[,c(2,3,4)]
head(count_sub )

gene_length <- count[,c(1)]
head(gene_length)

d <- DGEList(counts=count_sub,genes = gene_length)
d <- calcNormFactors(d, method="TMM")

rpkm<-rpkm(d,gene_length)
rpkm<-as.data.frame(rpkm)

rpkm$Name<-row.names(rpkm)
rpkm <- rpkm %>% group_by(Name) %>% mutate(avg=mean(NumReads,NumReads.1,NumReads.2))
rpkm <- rpkm %>% filter(avg>0) %>% mutate(trans=log2(avg))

pdf("Expression_Isr7_final.pdf", height=5, width = 10)
ggplot(rpkm) + geom_density(aes(x=trans), color="black",  fill="orange", alpha=0.5) +
  geom_vline(xintercept=c(median(rpkm$trans),mean(rpkm$trans),quantile(rpkm$trans, 0.90),quantile(rpkm$trans, 0.95)), linetype='dashed', color="black") + theme_bw() + xlim(min(rpkm$trans),max(rpkm$trans)) +
  geom_point(aes(x=rpkm$trans[rpkm$Name=="BgISR7-10067"], y=0)) + ylab("density") + xlab("log2(rpkm)")
dev.off()
