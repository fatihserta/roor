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



#NMDS analysis of PA table

    library(vegan)
    library(ggplot2)

    nmds <- metaMDS(PA,distance = "jaccard",k = 3)
    
    nmds$stress
    stressplot(nmds)
    plot(nmds, display = c("sites"),type = "n")
    text(nmds, display = c("sites"),labels = row.names(PA),col="red")
    

#Fitness of environmental variables
    
    
    
    #vf<-envfit(nmds,ENV,perm = 999,choices=c(1,2))
    
    #significant variables
    vf<-envfit(nmds,ENV[,c(3:5,8,13,18,24,30,32:36,39:43)],perm = 999,choices=c(1,2))
    vf
    
#sample and region names 
       
    samp.name <- cbind(c("PKF-1-2",	"PKF-1-3",	"PKF-1-4",	"PKF-2-1",	"PKF-2-2",	"PKF-2-3",	
                         "PKF-2-4",	"PKF-3-1",	"PKF-3-2",	"PKF-3-3",	"PKF-3-4",	"PKF-4-1",	
                         "PKF-4-2",	"PKF-4-3",	"PKF-4-4",	"PKF-5-1",	"PKF-5-2",	"PKF-5-3",	
                         "PKF-5-4",	"PKF-6-1",	"PKF-6-2",	"PKF-6-3",	"PKF-6-4",	"PKF-7-1",	
                         "PKF-7-2",	"PKF-7-3",	"PKF-7-4",	"PKF-8-1",	"PKF-8-2",	"PKF-8-3",	
                         "PKF-8-4",	"YP-1",	"YP-2",	"YP-3",	"YP-4",	"YP-5",	"YP-6",	"VR-1",	
                         "VR-2",	"VR-3",	"VR-4",	"VR-5",	"SS-1-1",	"SS-1-2",	"SS-1-3",	"SS-1-4",	
                         "SS-1-5",	"SS-2-1",	"SS-2-2",	"SS-2-3",	"SS-2-4",	"SS-2-5",	"SP-1-2",	
                         "SP-1-3",	"SP-1-4",	"SP-2-1",	"SP-2-2",	"SP-2-3",	"SP-2-4",	"SP-3-1",	
                         "SP-3-2",	"SP-3-3",	"SP-3-4",	"OB-1-1",	"OB-1-2",	"OB-1-3",	"OB-1-4",	
                         "OB-2-1",	"OB-2-2",	"OB-2-3",	"OB-2-4",	"OB-2-5",	"OB-2-6",	"OB-3-1",	
                         "OB-3-2",	"OB-3-3",	"OB-3-4",	"OB-3-5"),ENV$seep)
    
    reg.name <-c(rep("PKF-S",15),rep("PKF-NS",12),rep("PKF-S",4),
                 rep("YP+VR",11),rep("SP+SS",21),rep("OB",15))
    
    
    
# nmds score table    
    scrs <- as.data.frame(scores(nmds, display = "sites"))
    scrs <- cbind(scrs, Group = reg.name)
    
# vector table
    spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
    spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
    
    
# covariance ellipse function
    
    veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))}

# coveriance ellipse plot table
    
    
    ord<-ordiellipse(nmds, as.factor(reg.name), display = "sites",
                     kind = "sd", conf = .95, label = F)
     
    df_ell <- data.frame()
    for(g in levels(scrs$Group)){
      df_ell <- rbind(df_ell,
                      cbind(as.data.frame(with(scrs[scrs$Group==g,],
                                               veganCovEllipse(ord[[g]]$cov,ord[[g]]$center))),
                            Group=g))}
 
    
    
    
    # plot based on nmds values. Point colors are given based on regions. 
    # Seep vs. non-seep are given as text color. Black: non-seep and red: seep. 
    
    pcol<- c("navyblue","#3B9AB2", "#9EBE91", "#E4B80E", "#F21A00")    
       
    ggplot(scrs, aes(x = NMDS1, y = NMDS2))+
      
      #geom_point(data=scrs2, aes(x = NMDS1, y = NMDS2),shape=16,color="lightgrey",size=0.2)+
      geom_point(aes( fill= factor(Group)),alpha=1,shape=21,size=4) +
      geom_text_repel(label=samp.name[,1],cex=3,hjust = 0,nudge_x = 0.008,nudge_y = 0.008,
                      color=factor(samp.name[,2]))+
      guides(alpha=FALSE,fill=guide_legend(title = element_blank()),
             color=guide_legend(title = element_blank()))+
      scale_fill_manual(values=pcol)+
      
      theme(legend.justification=c(0,0), legend.position=c(0.02,0),
            legend.text=element_text(size=10,face="italic"),legend.key=element_blank())+
      coord_fixed()+  ## need aspect ratio of 1!
      #vectors
      geom_segment(data = spp.scrs,
                   aes(x = 0, xend = NMDS1*0.8, y = 0, yend = NMDS2*0.8),
                   arrow = arrow(length = unit(0.2, "cm")), colour = "grey")+
      geom_text(data = spp.scrs, 
                aes(x = NMDS1*0.85, y = NMDS2*0.85, label = Species), 
                size = 4,color="darkgrey")+
      
      theme(panel.border = element_rect(fill=NA,color="black", linetype="solid"),
            panel.background = element_blank())+
      theme(axis.text.x = element_text(colour="black",size=12),
            axis.text.y = element_text(colour="black",size=12),  
            axis.title.x = element_text(colour="black",size=12 ),
            axis.title.y = element_text(colour="black",size=12))+
      theme(panel.grid.major= element_blank(), 
            panel.grid.minor = element_blank())+
      theme(plot.title = element_text(lineheight=1.8, face="bold", color = "red"))+
      geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,color=factor(Group)),  linetype=1)+
      scale_color_manual(values=pcol,name = "St.Dev.\nEllipse")+
      #guides(alpha=FALSE,fill=guide_legend(title="St.Dev.\nEllipse"))+
      #xlim(-1,0.6)+
      #ylim(-1,0.7)+
      annotate("text", x = 0.38, y = -0.75, label = "Stress = 0.08", fontface="italic")
    