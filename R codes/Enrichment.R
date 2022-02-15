library(ReactomePA)#
library(clusterProfiler)#
library(ggplot2)#
library(openxlsx)
library(viridis)
library(stringr)


mainDir <-"~/folder/path/input"
outDir <- "~/folder/path/output"
dir.create(file.path(outDir))

#IgG
cond <- ("igg")
dir.create(file.path(outDir, cond))
#Run Enrichment
geneList <- "input.data"

edo <- enrichPathway(geneList$IgG_entrez)
egoCC <- enrichGO(geneList$IgG_entrez, ont="CC",  OrgDb = "org.Hs.eg.db")
egoBP <- enrichGO(geneList$IgG_entrez, ont="BP",  OrgDb = "org.Hs.eg.db")
egoMF <- enrichGO(geneList$IgG_entrez, ont="MF",  OrgDb = "org.Hs.eg.db")

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')
egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')

###############################################################
## PLOTS (Pahtways)

#Dot plot of enriched terms. ORA
png(file = file.path(outDir, cond, "igg_ORA.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(edo, 
        showCategory=10, 
        x = "GeneRatio",
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#fb6a4a", "#b10026")) + 
  ggtitle("IgG - dotplot for ORA") +
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 29))
dev.off()

###############################################################
## PLOTS (CC)

#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "igg_CC.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(egoCC, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#fb6a4a", "#b10026")) + 
  ggtitle("IgG - dotplot for CC") +
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 18))
dev.off()

###############################################################
## PLOTS (BP)
#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "igg_BP.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(egoBP, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#fb6a4a", "#b10026")) + 
  ggtitle("IgG - dotplot for BP") +
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 20))
dev.off()

###############################################################
## PLOTS (MF)

#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "igg_MF.png"), bg = "transparent", width = 2000, height = 900, units = "px", res = 300)
barplot(egoMF, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#fb6a4a", "#b10026")) + 
  ggtitle("IgG - dotplot for MF") +
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

###############################################################
#Export tables

tableORA <- edox@result
tableBP <- egoxBP@result
tableCC <- egoxCC@result
tableMF <- egoxMF@result

tabelaORA <- subset(tableORA, p.adjust < 0.05)
tabelaBP <- subset(tableBP, p.adjust < 0.05)
tabelaCC <- subset(tableCC, p.adjust < 0.05)
tabelaMF <- subset(tableMF, p.adjust < 0.05)

write.xlsx(as.data.frame(tableORA), file.path(outDir, cond, "igg_tabelaORA.xlsx"))
write.xlsx(as.data.frame(tableBP), file.path(outDir, cond, "igg_tabelaBP.xlsx"))
write.xlsx(as.data.frame(tableCC), file.path(outDir, cond, "igg_tabelaCC.xlsx"))
write.xlsx(as.data.frame(tableMF), file.path(outDir, cond, "igg_tabelaMF.xlsx"))

write.xlsx(as.data.frame(tabelaORA), file.path(outDir, cond, "igg_tabelaORA_sig.xlsx"))
write.xlsx(as.data.frame(tabelaBP), file.path(outDir, cond, "igg_tabelaBP_sig.xlsx"))
write.xlsx(as.data.frame(tabelaCC), file.path(outDir, cond, "igg_tabelaCC_sig.xlsx"))
write.xlsx(as.data.frame(tabelaMF), file.path(outDir, cond, "igg_tabelaMF_sig.xlsx"))


###################################################################################################################

#IgA

cond <- ("iga")
dir.create(file.path(outDir, cond))
#Run Enrichment
geneList <- read.xlsx(file.path(mainDir, "iga.xlsx"), colNames = T)

edo <- enrichPathway(geneList$IgA_entrez)
egoCC <- enrichGO(geneList$IgA_entrez, ont="CC",  OrgDb = "org.Hs.eg.db")
egoBP <- enrichGO(geneList$IgA_entrez, ont="BP",  OrgDb = "org.Hs.eg.db")
egoMF <- enrichGO(geneList$IgA_entrez, ont="MF",  OrgDb = "org.Hs.eg.db")

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')
egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')

###############################################################
## PLOTS (Pahtways)

#Dot plot of enriched terms. ORA
png(file = file.path(outDir, cond, "iga_ORA.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(edo, 
        showCategory=10, 
        x = "GeneRatio",
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#4eb3d3", "#0c2c84")) + 
  ggtitle("IgA - dotplot for ORA") +
  labs(x = "GeneRatio") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))
dev.off()

###############################################################
## PLOTS (CC)

#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "iga_CC.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(egoCC, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#4eb3d3", "#0c2c84")) +  
  ggtitle("IgA - dotplot for CC") +
  labs(x = "GeneRatio") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

###############################################################
## PLOTS (BP)
#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "iga_BP.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(egoBP, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#4eb3d3", "#0c2c84")) +  
  ggtitle("IgA - dotplot for BP") + 
  labs(x = "GeneRatio") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))

dev.off()

###############################################################
## PLOTS (MF)

#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "iga_MF.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(egoMF, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#4eb3d3", "#0c2c84")) +  
  ggtitle("IgA - dotplot for MF") +
  labs(x = "GeneRatio") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

###############################################################
#Export tables

tableORA <- edox@result
tableBP <- egoxBP@result
tableCC <- egoxCC@result
tableMF <- egoxMF@result

tabelaORA <- subset(tableORA, p.adjust < 0.05)
tabelaBP <- subset(tableBP, p.adjust < 0.05)
tabelaCC <- subset(tableCC, p.adjust < 0.05)
tabelaMF <- subset(tableMF, p.adjust < 0.05)

write.xlsx(as.data.frame(tableORA), file.path(outDir, cond, "iga_tabelaORA.xlsx"))
write.xlsx(as.data.frame(tableBP), file.path(outDir, cond, "iga_tabelaBP.xlsx"))
write.xlsx(as.data.frame(tableCC), file.path(outDir, cond, "iga_tabelaCC.xlsx"))
write.xlsx(as.data.frame(tableMF), file.path(outDir, cond, "iga_tabelaMF.xlsx"))

write.xlsx(as.data.frame(tabelaORA), file.path(outDir, cond, "iga_tabelaORA_sig.xlsx"))
write.xlsx(as.data.frame(tabelaBP), file.path(outDir, cond, "iga_tabelaBP_sig.xlsx"))
write.xlsx(as.data.frame(tabelaCC), file.path(outDir, cond, "iga_tabelaCC_sig.xlsx"))
write.xlsx(as.data.frame(tabelaMF), file.path(outDir, cond, "iga_tabelaMF_sig.xlsx"))



###################################################################################################################


#IgA_IgG

cond <- ("iga_igg")
dir.create(file.path(outDir, cond))
#Run Enrichment
geneList <- read.xlsx(file.path(mainDir, "igs.xlsx"), colNames = T)

edo <- enrichPathway(geneList$IgA_entrez)
egoCC <- enrichGO(geneList$IgA_entrez, ont="CC",  OrgDb = "org.Hs.eg.db")
egoBP <- enrichGO(geneList$IgA_entrez, ont="BP",  OrgDb = "org.Hs.eg.db")
egoMF <- enrichGO(geneList$IgA_entrez, ont="MF",  OrgDb = "org.Hs.eg.db")

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')
egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')

###############################################################
## PLOTS (Pahtways)

#Dot plot of enriched terms. ORA
png(file = file.path(outDir, cond, "iga_igg_ORA.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(edo, 
        showCategory=10, 
        x = "GeneRatio",
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#b10026" ,"#fb6a4a", "#4eb3d3", "#0c2c84")) + 
  ggtitle("IgA_IgG - dotplot for ORA") +
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

###############################################################
## PLOTS (CC)

#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "iga_igg_CC.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(egoCC, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#b10026" ,"#fb6a4a", "#4eb3d3", "#0c2c84")) + 
  ggtitle("IgA_IgG - dotplot for CC") +
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

###############################################################
## PLOTS (BP)
#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "iga_igg_BP.png"), bg = "transparent", width = 2000, height = 1800, units = "px", res = 300)
barplot(egoBP, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#b10026" ,"#fb6a4a", "#4eb3d3", "#0c2c84")) + 
  ggtitle("IgA_IgG - dotplot for BP") + 
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

###############################################################
## PLOTS (MF)

#Dot plot of enriched terms.
png(file = file.path(outDir, cond, "iga_igg_MF.png"), bg = "transparent", width = 2000, height = 1000, units = "px", res = 300)
barplot(egoMF, 
        showCategory=10, 
        x = "GeneRatio", 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#b10026" ,"#fb6a4a", "#4eb3d3", "#0c2c84")) + 
  ggtitle("IgA_IgG - dotplot for MF") + 
  labs(x = "GeneRatio")+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30))
dev.off()

###############################################################
#Export tables

tableORA <- edox@result
tableBP <- egoxBP@result
tableCC <- egoxCC@result
tableMF <- egoxMF@result

tabelaORA <- subset(tableORA, p.adjust < 0.05)
tabelaBP <- subset(tableBP, p.adjust < 0.05)
tabelaCC <- subset(tableCC, p.adjust < 0.05)
tabelaMF <- subset(tableMF, p.adjust < 0.05)

write.xlsx(as.data.frame(tableORA), file.path(outDir, cond, "name.output"))
write.xlsx(as.data.frame(tableBP), file.path(outDir, cond, "name.output"))
write.xlsx(as.data.frame(tableCC), file.path(outDir, cond, "name.output"))
write.xlsx(as.data.frame(tableMF), file.path(outDir, cond, "name.output"))

write.xlsx(as.data.frame(tabelaORA), file.path(outDir, cond, "name.output"))
write.xlsx(as.data.frame(tabelaBP), file.path(outDir, cond, "name.output"))
write.xlsx(as.data.frame(tabelaCC), file.path(outDir, cond, "name.output"))
write.xlsx(as.data.frame(tabelaMF), file.path(outDir, cond, "name.output"))


barplot(egoBP, showCategory=10) + scale_fill_viridis()
barplot(egoBP, showCategory=10) + scale_fill_gradient2(low = "#a8ddb5", mid = "#4eb3d3", high = "#0c2c84")
barplot(egoBP, showCategory=10) + scale_fill_gradient2(low = "#a8ddb5", mid = "#fb6a4a", high = "#b10026")
barplot(egoBP, showCategory=10) + scale_fill_gradient(low = "#a8ddb5", high = "#fb6a4a")
barplot(egoBP, showCategory=10) + scale_fill_gradient(low = "#a8ddb5", high = "#b10026")
barplot(egoBP, showCategory=10) + scale_fill_gradient(low = "#fb6a4a", high = "#b10026")

#IgG
barplot(egoBP, 
        showCategory=10, 
        font.size = 12) + 
  scale_fill_gradientn(colours = c("#a8ddb5", "#fb6a4a", "#b10026")) + 
  ggtitle("IgG - dotplot for BP")
#IgA
barplot(egoBP, showCategory=10) + scale_fill_gradientn(colours = c("#a8ddb5", "#4eb3d3", "#0c2c84"))
#IgA_IgA
barplot(egoBP, showCategory=10) + scale_fill_gradientn(colours = c("#4eb3d3", "#0c2c84", "#b10026"))
barplot(egoBP, showCategory=10) + scale_fill_gradientn(colours = c("#b10026", "#4eb3d3", "#0c2c84"))
barplot(egoBP, showCategory=10) + scale_fill_gradientn(colours = c("#fb6a4a", "#0c2c84", "#4eb3d3"))
barplot(egoBP, showCategory=10) + scale_fill_gradientn(colours = c("#b10026" ,"#fb6a4a","#a8ddb5", "#4eb3d3", "#0c2c84"))
barplot(egoBP, showCategory=10) + scale_fill_gradientn(colours = c("#b10026" ,"#fb6a4a", "#4eb3d3", "#0c2c84"))

scale_x_discrete(labels = function(x) str_wrap(x, width = 30))







