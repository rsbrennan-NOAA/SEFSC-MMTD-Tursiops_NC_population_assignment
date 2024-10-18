#-------------------------------------------------------------------------------
# DAPC
library(adegenet)
library(dartR)
library(ggpubr)

dat <- read.structure(
  file= "./data/6479snps_run4_K4structpops_modified.str",
  onerowperind = FALSE,
  n.loc=6479 ,
  n.ind=136 ,
  col.lab = 0,
  col.pop = 2,
  NA.char = "-9",
  col.others=1,
  row.marknames=1
)

genl <- gi2gl(dat, parallel = FALSE, verbose = NULL)
genl@ind.names <- dat@other$X

# PCA

pca1 <- glPca(genl,center = T, scale = T, nf = 5)

png("figures/pca_eigen.png", h=4, w=4, units="in", res=300)
barplot(100*pca1$eig/sum(pca1$eig), col = heat.colors(50), main="PCA Eigenvalues") # retain first 5 axes, incremental decrease after 2
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
dev.off()

#proportion of explained variance by first four axes
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
a4<-pca1$eig[4]/sum(pca1$eig) # proportion of variation explained by 4th axis
pcvar <- data.frame(Axis = c(1:4), Proportion = c(a1,a2,a3,a4))
pcvar

pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$pop <- pop(genl)
pca1.scores$ind <- genl@ind.names

set.seed(89)
num_pops <- length(levels(factor(pca1.scores$pop)))

# plot PC 1 and 2

p1 <- ggplot(pca1.scores, aes(x=PC1, y=PC2, fill=pop))+
  geom_point(size=4, pch=21)+
  xlab(paste0("PC1 (",round(pcvar[1,2]*100,2),"%)") )+
  ylab(paste0("PC2 (",round(pcvar[2,2]*100,2),"%)")) + 
  theme_bw(base_size = 14) +
  scale_fill_manual(values=c("orange3", "red3", "darkgreen", "lawngreen"))
p1

ggsave("figures/pca.png", p1, h=4, w=6)

############ ----------------------------------------------------------------

# follow the k-1 recommendation for the number of PCAs. So k=4
grp <- find.clusters(genl, n.pca=3, max.n.clust=40)
# 4 is where it plateaus. maybe 3

table.value(table(pop(genl), grp$grp), col.lab=paste("inf", 1:6))


# Calculate optimal number of PCs 
dapcTemp <- dapc(genl, genl@pop, 
                 n.pca=60, n.da = 4)   

ascores <- optim.a.score(dapcTemp, smart = FALSE, n.sim = 50)
# 3 is best. T

library("poppr")
#pramx <- xvalDapc(tab(genl, NA.method = "mean"), pop(genl)) # says 10.
pramx <- xvalDapc(tab(genl, NA.method = "mean"), grp$grp) # says 10.

names(pramx)
pramx[-1]

scatter(pramx$DAPC, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

dapc2 <- dapc(genl, , 
              n.pca=3, n.da = 3)
scatter(dapc2)

###### --------------------------------------------------------------------------
# leave one out, assignment success

# leave out 30% of samples.
kept.id <- sample(1:nInd(genl), size=round(nInd(genl)*0.7), replace=F)

x <- genl[kept.id,]
x.sup <- genl[!1:nInd(genl) %in% kept.id,]

dapc4 <- dapc(x ,n.pca=3,n.da=3)
pred.sup <- predict.dapc(dapc4, newdata=x.sup )
names(pred.sup)

# assignments of indivs
pred.sup$assign
#coords of scores:
pred.sup$ind.scores
# posterior membership probs
round(pred.sup$posterior[1:5, 1:4],3)

# plot results:
col <- rainbow(length(levels(pop(x))))
col.points <- transp(col[as.integer(pop(x))],.2)
scatter(dapc4, col=col, bg="white", scree.da=0, pch="",
        cstar=0, clab=0, legend=TRUE)
par(xpd=TRUE)
points(dapc4$ind.coord[,1], dapc4$ind.coord[,2], pch=20,
       col=col.points, cex=5)
col.sup <- col[as.integer(pop(x.sup))]
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,2], pch=15,
       col=transp(col.sup,.9), cex=2)
add.scatter.eig(dapc4$eig,15,1,2, posi="bottomright", inset=.02)

mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))

table(pred.sup$assign, pop(x.sup))

# loop to get some estimate of accuracy
n=100
acc.out.4 <- rep(NA, n)
df.out.4 <- as.data.frame(matrix(nrow=n, ncol=4))
colnames(df.out.4) <- c("Orange_84", "Red_21", "DarkGreen_17", "LightGreen_14")

for(i in 1:n){
  # leave out ~20% of samples.
  kept.id <- sample(1:nInd(genl), size=round(nInd(genl)*0.8), replace=F)
  x <- genl[kept.id,]
  x.sup <- genl[!(1:nInd(genl)) %in% kept.id,]
  # run dapc
  dapc4 <- dapc(x , x@pop,n.pca=4,n.da=5)
  pred.sup <- predict.dapc(dapc4, newdata=x.sup )
  
  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out.4[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  
  # get accuracy by population
  #
  orange.tmp <- which(as.character(pop(x.sup)) == "Orange_84")
  if(length(orange.tmp) > 0){
    df.out.4$Orange_84[i] <- mean(as.character(pred.sup$assign[orange.tmp])==as.character(pop(x.sup[orange.tmp])))
  }
  red.tmp <- which(as.character(pop(x.sup)) == "Red_21")
  if(length(red.tmp) > 0){
    df.out.4$Red_21[i] <- mean(as.character(pred.sup$assign[red.tmp])==as.character(pop(x.sup[red.tmp])))
  }
  darkgreen.tmp <- which(as.character(pop(x.sup)) == "DarkGreen_17")
  if(length(darkgreen.tmp) > 0){
    df.out.4$DarkGreen_17[i] <- mean(as.character(pred.sup$assign[darkgreen.tmp])==as.character(pop(x.sup[darkgreen.tmp])))
  }
  lightgreen.tmp <- which(as.character(pop(x.sup)) == "LightGreen_14")
  if(length(lightgreen.tmp) > 0){
    df.out.4$LightGreen_14[i] <- mean(as.character(pred.sup$assign[lightgreen.tmp])==as.character(pop(x.sup[lightgreen.tmp])))
  }
  print(i)
}

write.csv(file="analysis/dapc_accuracy_prop.csv", df.out.4, row.names=F)

library(tidyr)

df_long <- df.out.4 %>%
  pivot_longer(
    cols = everything(),
    names_to = "population",
    values_to = "accuracy"
  )


p1 <- ggplot(df_long, aes(x=population, y=accuracy, fill=population))+
  geom_boxplot()+
  ylab("% correct assignments") +
  xlab("population") + 
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none") +
  ylim(0,1) +
  scale_fill_manual(values=c("darkgreen", "lawngreen", "orange3", "red3")) +
  ggtitle("DAPC: all individuals")

ggsave("figures/dapc_allindivs_accuracy.png",p1,
       h=4, w=5)



#-------------------------------------------------------------------------------
# downsample to 14 indivs, leave out 4 indivs

# loop to get some estimate of accuracy
n=100
acc.out.4 <- rep(NA, n)
df.out.fixed <- as.data.frame(matrix(nrow=n, ncol=4))
colnames(df.out.fixed) <- c("Orange_84", "Red_21", "DarkGreen_17", "LightGreen_14")

orange <- which(genl@pop == "Orange_84")
Lgreen <- which(genl@pop == "LightGreen_14")
Dgreen <- which(genl@pop == "DarkGreen_17")
red   <- which(genl@pop == "Red_21")
sampN <- 10

for(i in 1:n){

  orange.in <- sample(orange, sampN, replace=F)
  Lgreen.in <- sample(Lgreen, sampN, replace=F)
  Dgreen.in <- sample(Dgreen, sampN, replace=F)
  red.in    <- sample(red, sampN, replace=F)
  kept.id <- c(orange.in,Lgreen.in,Dgreen.in,red.in )
  x.sup <- genl[!(1:nInd(genl)) %in% kept.id,]
  x <- genl[kept.id,]
  
  # run dapc
  dapc4 <- dapc(x , x@pop,n.pca=4,n.da=5)
  pred.sup <- predict.dapc(dapc4, newdata=x.sup )
  
  # calculate accuracy of assignments. i.e., is anything coming from its actual population
  acc.out.4[i] <- mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))
  
  # get accuracy by population
  #
  orange.tmp <- which(as.character(pop(x.sup)) == "Orange_84")
  if(length(orange.tmp) > 0){
    df.out.fixed$Orange_84[i] <- mean(as.character(pred.sup$assign[orange.tmp])==as.character(pop(x.sup[orange.tmp])))
  }
  red.tmp <- which(as.character(pop(x.sup)) == "Red_21")
  if(length(red.tmp) > 0){
    df.out.fixed$Red_21[i] <- mean(as.character(pred.sup$assign[red.tmp])==as.character(pop(x.sup[red.tmp])))
  }
  darkgreen.tmp <- which(as.character(pop(x.sup)) == "DarkGreen_17")
  if(length(darkgreen.tmp) > 0){
    df.out.fixed$DarkGreen_17[i] <- mean(as.character(pred.sup$assign[darkgreen.tmp])==as.character(pop(x.sup[darkgreen.tmp])))
  }
  lightgreen.tmp <- which(as.character(pop(x.sup)) == "LightGreen_14")
  if(length(lightgreen.tmp) > 0){
    df.out.fixed$LightGreen_14[i] <- mean(as.character(pred.sup$assign[lightgreen.tmp])==as.character(pop(x.sup[lightgreen.tmp])))
  }
  print(i)
}

write.csv(file="analysis/dapc_accuracy_fixed.csv", df.out.fixed, row.names=F)

df_long <- df.out.fixed %>%
  pivot_longer(
    cols = everything(),
    names_to = "population",
    values_to = "accuracy"
  )


p1 <- ggplot(df_long, aes(x=population, y=accuracy, fill=population))+
  geom_boxplot()+
  ylab("% correct assignments") +
  xlab("population") + 
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none") +
  ylim(0,1) +
  scale_fill_manual(values=c("darkgreen", "lawngreen", "orange3", "red3")) +
  ggtitle("DAPC: 10 individual training")

ggsave("figures/dapc_fixedIndiv_accuracy.png",p1,
       h=4, w=5)

