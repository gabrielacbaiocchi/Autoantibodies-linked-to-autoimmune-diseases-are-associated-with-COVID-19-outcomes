# libraries:
rm(list=ls())

my_packages <- c("readxl","MSFA", "ggplot2", "ggpubr")
lapply(my_packages, library, character.only = TRUE)
 

# Data: Data only contains the columns required for the BMSFA.
path <- "cData_all_abs_new_gname.xlsx"
data_antibodies <- read_excel(path = path, sheet = "All")
group.name <- unique(data_antibodies$Group)
group.name
names(data_antibodies)

# Split data into four groups and get together. 
# Remove the first two columns:
data.control <- data_antibodies[data_antibodies$Group == "Control", -c(1:2)] 
data.mild <- data_antibodies[data_antibodies$Group == "Mild", -c(1:2)]
data.moderate <- data_antibodies[data_antibodies$Group == "Moderate", -c(1:2)]
data.severe <- data_antibodies[data_antibodies$Group == "Severe", -c(1:2)]

X_s <- list(data.control, data.mild, data.moderate, data.severe)

# centered data:
cX_s <- lapply(X_s, scale, center = TRUE, scale = FALSE) # X_s in sp_msfa will be automatically scaled.


# Fitting the BMSFA:
k = 5; j_s <- c(5, 5, 5, 5) 

fit.bmsfa <- sp_msfa(X_s = cX_s, k = k, j_s = j_s, outputlevel = 1)



nchain <- dim(fit.bmsfa$Phi)[3]
p <- dim(fit.bmsfa$Phi)[1]

Sig_phi <- array(NA, dim=c(p, p, nchain))
Sig_lcontrol <- array(NA, dim=c(p, p, nchain))
Sig_lmild <- array(NA, dim=c(p, p, nchain))
Sig_lmoderate <- array(NA, dim=c(p, p, nchain))
Sig_lsevere <- array(NA, dim=c(p, p, nchain))
 

# Big matrix requires memory here:
for(i in 1:nchain){
  Sig_phi[,,i] = crossprod(t(fit.bmsfa$Phi[,,i]))
  Sig_lcontrol[,,i] = crossprod(t(fit.bmsfa$Lambda[[1]][,,i]))
  Sig_lmild[,,i] = crossprod(t(fit.bmsfa$Lambda[[2]][,,i]))
  Sig_lmoderate[,,i] = crossprod(t(fit.bmsfa$Lambda[[3]][,,i]))
  Sig_lsevere[,,i] = crossprod(t(fit.bmsfa$Lambda[[4]][,,i]))
}
	

# Estimated cov matrices:
Sig_phiE = apply(Sig_phi , c(1,2), mean)
Sig_LcontrolE = apply(Sig_lcontrol , c(1,2), mean)
Sig_LmildE = apply(Sig_lmild , c(1,2), mean)
Sig_LmoderateE = apply(Sig_lmoderate , c(1,2), mean)
Sig_LsevereE = apply(Sig_lsevere , c(1,2), mean)


# Selection for the number of common factor:
fun_eigen <- eigen(Sig_phiE)
val_eigen <- fun_eigen$values
prop_var <- val_eigen/sum(val_eigen)
choose_K <- length(which(prop_var > 0.05))

# Selection for the number of study-specific factors:
fun_eigenControl <- eigen(Sig_LcontrolE)
val_eigenControl <- fun_eigenControl$values
prop_varControl <- val_eigenControl/sum(val_eigenControl)
choose_JControl <- length(which(prop_varControl > 0.05))

fun_eigenMild <- eigen(Sig_LmildE)
val_eigenMild <- fun_eigenMild$values
prop_varMild <- val_eigenMild/sum(val_eigenMild)
choose_JMild <- length(which(prop_varMild > 0.05))

fun_eigenModerate <- eigen(Sig_LmoderateE)
val_eigenModerate <- fun_eigenModerate$values
prop_varModerate <- val_eigenModerate/sum(val_eigenModerate)
choose_JModerate <- length(which(prop_varModerate > 0.05))

fun_eigenSevere <- eigen(Sig_LsevereE)
val_eigenSevere <- fun_eigenSevere$values
prop_varSevere <- val_eigenSevere/sum(val_eigenSevere)
choose_JSevere <- length(which(prop_varSevere > 0.05))
 

##########################################
# compute factor loading.
p_name <- colnames(cX_s[[1]])
 
k = choose_K
fun_eigen <- eigen(Sig_phiE)
val_matrix <- diag(sqrt(fun_eigen$values))
load <- fun_eigen$vec %*% val_matrix
loadK <- load[,1:k]
rownames(loadK) <- p_name

j_co = choose_JControl
j_mild = choose_JMild
j_moderate = choose_JModerate
j_severe = choose_JSevere


fun_eigenCO <- eigen(Sig_LcontrolE)
val_eigenCO <- diag(sqrt(fun_eigenCO$values))
l_co <- fun_eigenCO$vec %*% val_eigenCO
l_co <- l_co[,1:j_co]

fun_eigenMILD <- eigen(Sig_LmildE)
val_eigenMILD <- diag(sqrt(fun_eigenMILD$values))
l_mild <- fun_eigenMILD$vec %*% val_eigenMILD
l_mild <- l_mild[,1:j_mild]

fun_eigenMODERATE <- eigen(Sig_LmoderateE)
val_eigenMODERATE <- diag(sqrt(fun_eigenMODERATE$values))
l_moderate <- fun_eigenMODERATE$vec %*% val_eigenMODERATE
l_moderate <- l_moderate[,1:j_moderate]


fun_eigenSEVERE <- eigen(Sig_LsevereE)
val_eigenSEVERE <- diag(sqrt(fun_eigenSEVERE$values))
l_severe <- fun_eigenSEVERE$vec %*% val_eigenSEVERE
l_severe <- l_severe[,1:j_severe]


# Loadings matrices:
# loadK; l_co; l_mild; l_moderate; l_severe;


# PLOTS:
k <- dim(loadK)[2]
p <- dim(loadK)[1]

phi_com <- data.frame(Row = rep(p_name[1:p], times = k),  
                      Col = rep(x = c(" 1"," 2"), each = p), Y=matrix(c(loadK), p*k, 1))

library(ggplot2)
pD1 <- ggplot(phi_com, aes(Col , Row)) +ggtitle('Common')
pD1 <- pD1 + scale_y_reverse()
pD1 <- pD1 +   theme( legend.position = "none", 
                      plot.title = element_text(hjust = 0.5, vjust=2.12, size=11),
                      axis.text.y=element_text(size=2.5), axis.text.x=element_text(size=10),
                      axis.ticks=element_blank(), axis.title = element_text(size=10), 
                      panel.background = element_blank(), plot.margin=unit(c(1,0,1,-2.75),
                                                                           "lines"))
pD1 <- pD1 + scale_x_discrete("",labels=c(" 1"= expression(paste(phi)[1])," 2"= expression(paste(phi)[2])))
pD1 <- pD1 + scale_y_discrete("", limits=p_name[p:1])
pD1 <- pD1 + geom_tile(aes(fill=Y), colour="white") 
pD1 <- pD1 + scale_fill_gradient2(low = "#a8ddb5",
                                  mid = "white",
                                  high = "#878787", limits=c(-1, 1), '')


j_con=dim(l_co)[2]
l_coP <- data.frame(Row = rep(p_name[1:p], times = j_con),  
                    Col = rep(x = c(" 1"," 2"," 3"," 4"), each = p), Y=matrix(c(l_co), p*j_con, 1))

pD_co <- ggplot(l_coP, aes(Col , Row)) +ggtitle('Control')
pD_co <- pD_co + scale_y_reverse()
pD_co <- pD_co +   theme( legend.position = "none", plot.title = element_text(hjust = 0.5, vjust=2.12, size=11),
                          axis.text.y=element_blank(),axis.text.x=element_text(size=10),
                          axis.ticks=element_blank(), axis.title = element_text(size=10), 
                          panel.background = element_blank(), plot.margin=unit(c(1,0,1,-0.2), "lines"))
pD_co <- pD_co + scale_x_discrete("",labels=c(" 1"= expression(paste(lambda)[11])," 2"= expression(paste(lambda)[12]),
                                              " 3"=  expression(paste(lambda)[13])," 4"=  expression(paste(lambda)[14])))
pD_co <- pD_co + scale_y_discrete("", limits=p_name[p:1])
pD_co <- pD_co + geom_tile(aes(fill=Y), colour="white") 
pD_co <- pD_co + scale_fill_gradient2(low = "#a8ddb5",
                                      mid = "white",
                                      high = "#878787", limits=c(-1, 1), '')

j_mild=dim(l_mild)[2]
l_mildP <- data.frame(Row = rep(p_name[1:p], times = j_mild), 
                      Col = rep(x = c(" 1"," 2"," 3", " 4", " 5"), each = p), Y=matrix(c(l_mild), p*j_mild, 1))
pD_mild <- ggplot(l_mildP, aes(Col , Row)) +ggtitle('Mild')
pD_mild <- pD_mild + scale_y_reverse()
pD_mild <- pD_mild +   theme( legend.position = "none", plot.title = element_text(hjust = 0.5, vjust=2.12, size=11),
                              axis.text.y=element_blank(),axis.text.x=element_text(size=10),
                              axis.ticks=element_blank(), axis.title = element_text(size=10), 
                              panel.background = element_blank(), plot.margin=unit(c(1,0,1,-0.2), "lines"))
pD_mild <- pD_mild + scale_x_discrete("",labels=c(" 1"= expression(paste(lambda)[21])," 2"= expression(paste(lambda)[22]),
                                                  " 3"=  expression(paste(lambda)[23]), " 4"=  expression(paste(lambda)[24]), 
                                                  " 5"=  expression(paste(lambda)[25])))
pD_mild <- pD_mild + scale_y_discrete("", limits=p_name[p:1])
pD_mild <- pD_mild + geom_tile(aes(fill=Y), colour="white") 
pD_mild <- pD_mild+ scale_fill_gradient2(low = "#a8ddb5",
                                         mid = "white",
                                         high = "#878787", limits=c(-1, 1), '')


j_moderate = dim(l_moderate)[2]
l_moderateP <- data.frame(Row = rep(p_name[1:p], times = j_moderate), 
                          Col = rep(x = c(" 1"," 2"), each = p), Y=matrix(c(l_moderate), p*j_moderate, 1))
pD_moderate <- ggplot(l_moderateP, aes(Col , Row)) +ggtitle('Moderate')
pD_moderate <- pD_moderate + scale_y_reverse()
pD_moderate <- pD_moderate  +   theme( legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 2.12, size=11),axis.text.y=element_blank(),axis.text.x=element_text(size=10),axis.ticks=element_blank(), axis.title = element_text(size=10), panel.background = element_blank(), plot.margin=unit(c(1,0,1,-0.2), "lines"))
pD_moderate <- pD_moderate  + scale_x_discrete("",labels=c(" 1"= expression(paste(lambda)[31]),
                                                           " 2"= expression(paste(lambda)[32])))
pD_moderate <- pD_moderate  + scale_y_discrete("", limits=p_name[p:1])
pD_moderate <- pD_moderate + geom_tile(aes(fill=Y), colour="white") 
pD_moderate <- pD_moderate  + scale_fill_gradient2(low = "#a8ddb5",
                                                   mid = "white",
                                                   high = "#878787", limits=c(-1, 1), '')

j_severe = dim(l_severe)[2]
l_severeP <- data.frame(Row = rep(p_name[1:p], times = j_severe), 
                        Col = rep(x = c(" 1"," 2"), each = p), Y=matrix(c(l_severe), p*j_severe, 1))
pD_severe <- ggplot(l_severeP, aes(Col , Row)) +ggtitle('Severe')
pD_severe <- pD_severe + scale_y_reverse()
pD_severe <- pD_severe +   theme( legend.position = "none", plot.title = element_text(hjust = 0.5, vjust=2.12, size=11),axis.text.y=element_blank(),axis.text.x=element_text(size=10),axis.ticks=element_blank(), axis.title = element_text(size=10), panel.background = element_blank(), plot.margin=unit(c(1,0,1,-0.2), "lines"))
pD_severe <- pD_severe + scale_x_discrete("",labels=c(" 1"= expression(paste(lambda)[41]),
                                                      " 2"= expression(paste(lambda)[42])))
pD_severe <- pD_severe + scale_y_discrete("", limits=p_name[p:1])
pD_severe <- pD_severe + geom_tile(aes(fill=Y), colour="white") 
pD_severe <- pD_severe + scale_fill_gradient2(low = "#a8ddb5",
                                              mid = "white",
                                              high = "#878787", limits=c(-1, 1), '')

library(ggpubr)
pdf("Loadings_autoantibodies_corrected.pdf")
ggarrange(pD1, pD_co, pD_mild, pD_moderate, pD_severe, ncol=5, nrow=1, common.legend = T, legend = 'right')+
  theme(plot.margin = margin(-0.2,0.0,-0.75, 1.2, "cm"))  #margin. (t,r,b,l)
dev.off()



