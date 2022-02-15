library(factoextra)
library(ggExtra)
library(ggplot2)
library(scales)
library(reshape2)
library(viridis)

PCAf<- "input.data"
groups <- as.factor(PCAf$Group)

PCAf[1] <- NULL  

res.pca <- prcomp(PCAf, scale = TRUE)

#Biplot
p <- fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             geom = c("point"),
             palette = c("#939694", "#a8ddb5", "#4eb3d3", "#0c2c84"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
) + theme(legend.position = "bottom") + ggtitle("PCA") + 
theme(text = element_text(size = 25), axis.text = element_text(size = 18)) 
print(p)

#Adding histogram
p1 <- ggExtra::ggMarginal(p, type = "histogram", groupFill = TRUE, groupColour = TRUE)
print(p1)

#Biplot with vectors
p <- fviz_pca_biplot(res.pca,
                  col.ind = groups, # color by groups,
                  geom = c("point"),
                  col.var = "black",
                  alpha.var = 0.05,
                  palette = c("#939694", "#a8ddb5", "#4eb3d3", "#0c2c84"),
                  addEllipses = TRUE, # Concentration ellipses
                  ellipse.type = "confidence",
                  legend.title = "Groups: ",
                  repel = TRUE
) + theme(legend.position = "bottom") + ggtitle("PCA")+
theme(text = element_text(size = 25), axis.text = element_text(size = 20)) 
print(p)

#Adding histogram
p1 <- ggExtra::ggMarginal(p, type = "density", groupFill = TRUE, groupColour = TRUE)
print(p1)

#Heatmap
qplot(x=Var1, y=Var2, data=melt(cor(PCAf)), geom="tile",
                 fill=value) + #scale_fill_gradient2(
                   #low = "yellow",
                   #mid = "white",
                   #high = "royalblue",
                   #space = "Lab",
                   #na.value = "grey50",
                   #guide = "colourbar",
                   #aesthetics = "fill",
                   #limits=c(-1,1)
                 #) + 
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(colour = "black"))+ scale_fill_viridis(option="magma", limits=c(-1,1))


fviz_eig(res.pca)+
theme(text = element_text(size = 18), axis.text = element_text(size = 18))


fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#a6bddb", "#3690c0", "#034e7b"),
             alpha.var = 0.1,
             repel = TRUE     # Avoid text overlapping
)+
theme(text = element_text(size = 25), axis.text = element_text(size = 18))


fviz_pca_biplot(res.pca, repel = TRUE,
                alpha.var = 0.1,
                col.var = "#034e7b", # Variables color
                col.ind = "#696969"  # Individuals color
)+
theme(text = element_text(size = 18), axis.text = element_text(size = 18))

