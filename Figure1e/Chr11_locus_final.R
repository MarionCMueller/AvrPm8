####### Figure AvrPm8 Locus ######

library("dplyr")
library("ggplot2")
library("reshape2")
library(gggenomes)
library(thacklr)

#################################
### Final plot ##################
#################################



segChr11<-tibble(
  seq_id = c("ISR7_chr-11", "Bgt_chr-11"),
  length = c(400000,322122),
  start = c(2852000,2472626),
  end = c(3252000,2841868)
)

# a minimal gene track

g0chr11<-read.table("AvrPm8_gene_locus_comparison.txt",h=T)
g0chr11<-tibble(g0chr11)
g0chr11[1,2]<-g0chr11[1,2]-8000
g0chr11[2,3]<-g0chr11[2,3]-8000
g0chr11[3,2]<-g0chr11[3,2]+8000
g0chr11[4,2]<-g0chr11[4,2]-4000
g0chr11[4,3]<-g0chr11[4,3]+4000
g0chr11[5,3]<-g0chr11[5,3]-8000
g0chr11[6,3]<-g0chr11[6,3]+8000
g0chr11[7,3]<-g0chr11[7,3]+4000
g0chr11[7,2]<-g0chr11[7,2]-4000
g0chr11[8,2]<-g0chr11[8,2]-8000
g0chr11[9,2]<-g0chr11[9,2]-4000
g0chr11[9,3]<-g0chr11[9,3]+4000
g0chr11[10,3]<-g0chr11[10,3]+8000
g0chr11[11,2]<-g0chr11[11,2]-8000
g0chr11[12,2]<-g0chr11[12,2]-8000

g0chr11[13,2]<-g0chr11[13,2]-8000
g0chr11[14,3]<-g0chr11[14,3]-8000
g0chr11[15,2]<-g0chr11[15,2]+8000
g0chr11[16,2]<-g0chr11[16,2]-4000
g0chr11[16,3]<-g0chr11[16,3]+4000
g0chr11[17,3]<-g0chr11[17,3]-8000
g0chr11[18,3]<-g0chr11[18,3]+8000
g0chr11[19,3]<-g0chr11[19,3]+4000
g0chr11[19,2]<-g0chr11[19,2]-4000
g0chr11[20,2]<-g0chr11[20,2]-8000
g0chr11[21,2]<-g0chr11[21,2]-4000
g0chr11[21,3]<-g0chr11[21,3]+4000
g0chr11[22,3]<-g0chr11[22,3]+8000
g0chr11[23,2]<-g0chr11[23,2]-8000
g0chr11[24,2]<-g0chr11[24,2]-8000
is.factor(g0chr11$class)
is.factor(g0chr11$fam)
is.factor(g0chr11$strand)
g0chr11

link_data<-read.table("out.1coords",h=F)
link_data<-na.omit(link_data)
link_data
link_data_Chr11<- link_data  %>% select(V1,V2,V4,V3)
link_data_Chr11
link_data_chr11_avrpm8<-link_data %>% select(V12,V1,V2,V13,V3,V4)
names(link_data_chr11_avrpm8)<-c("seq_id","start","end","seq_id2","start2","end2")
link_data_chr11_avrpm8<-tibble(link_data_chr11_avrpm8)

borders_isr7<-tibble(
  seq_id = c("ISR7_chr-11","ISR7_chr-11"),
  start = c(2852000,3251000),
  end = c(2853000,3252000),
  name=c("2.85Mb","3.25Mb"))


borders_96224<-tibble(
  seq_id = c("Bgt_chr-11","Bgt_chr-11"),
  start = c(2472624,2840868),
  end = c(2473624,2841868),
  name=c("2.47Mb","2.84Mb"))

corr_factor<-3052007-2661079 

snpall<-read.table("AvrPm8_gene_locus_comparison.txt",h=T)
snps<-snpall[c(4,6,9,11),]
snps[,2]<-snps[,2]+500
snps[,3]<-snps[,3]-500
snps[2,2]<-snps[2,2]+4300
snps[2,3]<-snps[2,3]+4300
snps[4,2]<-snps[4,2]-4300
snps[4,3]<-snps[4,3]-4300
snps

snpsGwas<-tibble(
  seq_id = c("Bgt_chr-11"),
  start = c(3067561-300-corr_factor),
  end = c(3067561+300-corr_factor))
snpsGwas

pdf("LocusAvrPm8.pdf",w=20,h=5)
gggenomes(genes=g0chr11, seqs=segChr11, links=link_data_chr11_avrpm8,feats=list(borders_isr7,borders_96224,snps,snpsGwas))  +
  geom_feat(data=feats(borders_isr7),color="black",position="identity", size=5) +
  geom_feat_label(data=feats(borders_isr7),aes(label=name), nudge_y=0.1,check_overlap = TRUE,size=5,angle=90)  +
  geom_feat(data=feats(borders_96224),color="black",position="identity", size=5) +
  geom_feat_label(data=feats(borders_96224),aes(label=name), nudge_y=0.1,check_overlap = TRUE,size=5,angle=90,hjust=1.9) +
  geom_seq()  +      # draw contig/chromosome lines +
  geom_link(fill="grey", color="darkgrey") +
  geom_gene(aes(fill=fam),position="identity", size = 6, shape=c(6,10), intron_types=c()) +
  scale_fill_manual(values =c("#3B9AB2","grey","gold1"),labels=c('E003', 'E014', 'NonEff')) + theme(legend.position = "top",legend.title=element_blank())  +
  geom_gene_tag(aes(label=name), nudge_y=0.1, check_overlap = FALSE,size=4,angle=90) +
  geom_feat(data=feats(snps),color=c("red","black","black","black"),position="identity", size=10) +
  geom_feat(data=feats(snpsGwas),color=c("darkorange"),position="identity", size=10) +
  theme(legend.text=element_text(size=20))
dev.off()






