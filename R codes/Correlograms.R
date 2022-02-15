#CORRELOGRAM OF CIRCULES
#FROM: https://rpubs.com/bigcat/258548
install.packages("qgraph")
library(qgraph)
install.packages("openxlsx")
library(openxlsx)


IgA.all<- "input.data"

Control<-IgA.all[which(IgA.all[,1]=="Control"),1:42]
Mild<-IgA.all[which(IgA.all[,1]=="Mild"),1:42]
Severe<-IgA.all[which(IgA.all[,1]=="Severe"),1:42]
Oxygen<-IgA.all[which(IgA.all[,1]=="Oxygen"),1:42]

#Anosmia
Control.N.AN<-Control[which(Control$Anosmia=="Noanosmia"),1:42]
Mild.AN<-Mild[which(Mild$Anosmia=="Anosmia"),1:42]
Mild.N.AN<-Mild[which(Mild$Anosmia=="Noanosmia"),1:42]
Moderate.AN<-Severe[which(Severe$Anosmia=="Anosmia"),1:42]
Moderate.N.AN<-Severe[which(Severe$Anosmia=="Noanosmia"),1:42]
Severe.AN<-Severe[which(Oxygen$Anosmia=="Anosmia"),1:42]
Severe.N.AN<-Oxygen[which(Oxygen$Anosmia=="Noanosmia"),1:42]

#Removing symptoms
cor1<-Control.N.AN[,2:28]
cor2<-Mild.AN[,2:28]
cor3<-Mild.N.AN[,2:28]
cor4<-Moderate.AN[,2:28]
cor5<-Moderate.N.AN[,2:28]
cor6<-Severe.AN[,2:28]
cor7<-Severe.N.AN[,2:28]


write.xlsx(as.data.frame(IgG.all), "table.xlsx", row.names=TRUE)


par(mfrow=c(1,2))
View(cor1)


cormat=cor(cor1, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])
qgraph(cormat,shape="circle", posCol="#4eb3d3", negCol="#fb6a4a",
       layout="groups", vsize=8, label.cex=4, minimum=0.8, title = "Control No Anosmia")


cormat=cor(cor2, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])
qgraph(cormat,shape="circle", posCol="#4eb3d3", negCol="#fb6a4a",
       layout="groups", vsize=8, label.cex=4, minimum=0.8, title = "Mild Anosmia")

cormat=cor(cor3, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])
qgraph(cormat,shape="circle",  posCol="#4eb3d3", negCol="#fb6a4a",
       layout="groups", vsize=8, label.cex=4, minimum=0.8, title = "Mild No Anosmia")

cormat=cor(cor4, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])
qgraph(cormat,shape="circle",  posCol="#4eb3d3", negCol="#fb6a4a",
       layout="groups", vsize=8, label.cex=4, minimum=0.8, title = "Moderate Anosmia")

cormat=cor(cor5, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])
qgraph(cormat,shape="circle",  posCol="#4eb3d3", negCol="#fb6a4a",
       layout="groups", vsize=8, label.cex=4, minimum=0.8, title = "Moderate No Anosmia")

cormat=cor(cor6, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])
qgraph(cormat,shape="circle",  posCol="#4eb3d3", negCol="#fb6a4a",
       layout="groups", vsize=8, label.cex=4, minimum=0.8, title = "Severe Anosmia")

cormat=cor(cor7, method="spearman") # or cormat=cor(dat[,1:12],dat[,1:12])
qgraph(cormat,shape="circle",  posCol="#4eb3d3", negCol="#fb6a4a",
       layout="groups", vsize=8, label.cex=4, minimum=0.8, title = "Severe No Anosmia")


