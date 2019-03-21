#rm(list = ls())
#install.packages("repmis")

library(repmis)
source_data("https://github.com/fatihserta/roor/blob/master/combined.RData?raw=True")

PA[1:78,1:5] 
# Presence absence data (community table). Samples on rows molecular formulas on columns.
# if the forumula exist in that sample 1 if not 0.
# (eg: 10 11 1 4 1 means C10 H11 N O4 S)
ENV[1:78,1:5]
# Environmental variables and measured indexes for the samples




library(data.table)
library(vegan)
library(ggplot2)
library(FactoMineR)


#NMDS analysis of PA table

nmds <- metaMDS(PA,distance = "bray",k = 2)

nmds$stress
stressplot(nmds)
plot(nmds, display = c("sites"),type = "n")
text(nmds, display = c("sites"),labels = row.names(PA),col="red")

vf<-envfit(nmds,ENV,perm = 999,choices=c(1,2))
vf

scrs <- as.data.frame(scores(nmds, display = "sites"))
scrs <- cbind(scrs, Group = ENV$seep)
scrs <- cbind(scrs, Group = ENV$deep)
scrs <- cbind(scrs, Group = ENV$station)
scrs <- cbind(scrs, Group = ENV$mass)
scrs <- cbind(scrs, Group = ENV$metcat)

spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

ggplot(scrs, aes(x = NMDS1, y = NMDS2))+
  geom_label(aes( fill= factor(Group),alpha=0.3,label=rownames(PA)))  +
  guides(alpha=FALSE,fill=guide_legend(title="Water \nMasses"))+
  scale_fill_manual(values=c("darkgrey", "#E69F00", "#56B4E9","indianred","forestgreen"))+
  theme(legend.justification=c(0,0), legend.position=c(0.02,0),
        legend.text=element_text(size=10,face="italic"),legend.key=element_blank())+
  #coord_fixed()+  ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), colour = "darkgrey")+
  geom_text(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species), size = 3)+
  theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"))+
  theme(axis.text.x = element_text(colour="black",size=12),
        axis.text.y = element_text(colour="black",size=12),  
        axis.title.x = element_text(colour="black",size=12 ),
        axis.title.y = element_text(colour="black",size=12))+
  theme(panel.grid.major= element_blank(), 
        panel.grid.minor = element_blank())+
  theme(plot.title = element_text(lineheight=1.8, face="bold", color = "red"))
