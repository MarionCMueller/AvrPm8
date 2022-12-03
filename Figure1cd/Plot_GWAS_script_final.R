#################################################
######### GWAS Graph          ###################
#################################################

library(dplyr)
library(ggplot2)


#################################################################
######### Manhattan Plots AvrPm8 GWAS ###########################
#################################################################

################################
#### Pm8 OE #############
################################

gwas79Pm8OEStep<-read.csv("GWAS_79isolates/GAPIT.GLM.Pm8OE_2step.GWAS.Results.csv",h=T, sep=",",stringsAsFactors = FALSE)
head(gwas79Pm8OEStep)
nbmarkers3<-length(gwas79Pm8OEStep$P.value)
threshold3<--log10(0.05/nbmarkers3)
threshold3

head(gwas79Pm8OEStep)
gwas79Pm8OEStep <- gwas79Pm8OEStep %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-11" ,"11")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-10" ,"10")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-09" ,"9")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-08" ,"8")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-07" ,"7")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-06" ,"6")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-05" ,"5")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-04" ,"4")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-03" ,"3")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-02" ,"2")) %>%
  mutate(Chromosome=replace(Chromosome, Chromosome =="ISR7_chr-01" ,"1")) 
gwas79Pm8OEStep$Chromosome<-as.numeric(gwas79Pm8OEStep$Chromosome)
graph3 <- gwas79Pm8OEStep %>% group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  left_join(gwas79Pm8OEStep, ., by=c("Chromosome"="Chromosome")) %>%
  arrange(Chromosome, Position) %>%
  mutate(BPcum=Position+tot)
graph3

axisdf3 = graph3 %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
significandon3 <- graph3 %>% filter(-log10(P.value)>threshold3) # select significant SNPs
significandon3
ManPlot3<-ggplot() + geom_point(aes(x=graph3$BPcum, y=-log10(graph3$P.value),color=as.factor(graph3$Chromosome)), alpha=1, size=1.2) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(label = axisdf3$Chromosome, breaks= axisdf3$center ) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,9),breaks=seq(0,9,by=1.5)) +     
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
  ) +
  geom_hline(yintercept = threshold3, linetype=2, color="grey") + 
  geom_point(aes(x=significandon3$BPcum, y=-log10(significandon3$P.value)), color="black", fill=c("red","orange"), size=1.5, stroke=1, shape=21) +
  ylab("-log10(p.value)")

ManPlot3
pdf("GWASPm8OEGLM.pdf", width=10, height=2)
ManPlot3
dev.off()


################################
#### Pm8 OE Chr 11 #############
################################
graph3_chr11 <- graph3 %>% filter(Chromosome==11)
graph3_chr11

ManPlot3_chr11<-ggplot() + geom_point(aes(x=graph3_chr11$BPcum, y=-log10(graph3_chr11$P.value),color=as.factor(graph3_chr11$Chromosome)), alpha=1, size=1.1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(label = axisdf3$Chromosome, breaks= axisdf3$center ) +
  scale_y_continuous(expand = c(0, 0), limits=c(0,9),breaks=seq(0,9,by=1.5)) +     
  theme_classic() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x=element_blank(),
    axis.text.x=element_text(size=12),
    axis.text.y=element_text(size=12),
  ) +
  geom_hline(yintercept = threshold3, linetype=2, color="grey") + 
  geom_point(aes(x=significandon3$BPcum, y=-log10(significandon3$P.value)), color="black", fill=c("red","orange"), size=1.3, stroke=1, shape=21) +
  ylab("-log10(p.value)")

ManPlot3_chr11
pdf("GWASPm8OEGLM_chr11.pdf", width=3, height=2)
ManPlot3_chr11
dev.off()
