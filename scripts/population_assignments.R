

# summary of major approaches:
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13323
# https://repository.library.noaa.gov/view/noaa/55421/noaa_55421_DS1.pdf

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
  col.lab = NULL,
  col.pop = 2,
  NA.char = "-9",
  col.others=1,
  row.marknames=1
)




genl <- gi2gl(dat, parallel = FALSE, verbose = NULL)
genl@ind.names <- dat@other$X

# PCA

pca1 <- glPca(genl,center = T, scale = T, nf = 5)

#png("dapc_eigen.png", h=4, w=4, units="in", res=300)
barplot(100*pca1$eig/sum(pca1$eig), col = heat.colors(50), main="PCA Eigenvalues") # retain first 5 axes, incremental decrease after 2
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
#dev.off()

#proportion of explained variance by first three axes
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
pca1.p<-ggscatter(pca1.scores, x = "PC1", y = "PC2", color = "pop",
                  ellipse = T, ellipse.level = 0.95, size = 3,
                  xlab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"),
                  ylab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)")
)
pca1.p
#ggsave("dapc_PCA.png",pca1.p, h=5, w=5)
pca2.p<-ggscatter(pca1.scores, x = "PC2", y = "PC3", color = "pop",
                  ellipse = T, ellipse.level = 0.95,
                  xlab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  ylab = paste0("PC3 (",round(pcvar[3,2]*100,2),"%)"))
pca2.p



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

# leave out ~30% of samples.
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
  x.sup <- genl[!1:nInd(genl) %in% kept.id,]
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

boxplot(acc.out.4)

head(df.out.4)

boxplot(df.out.4)

png(filename = "./figures/dapc.accuracyPops.png", h=6, w=8, units="in", res=300)
boxplot(df.out.4, ylab="Proportion correct assignment")
dev.off()




















# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# assignPOP: https://alexkychen.github.io/assignPOP/analyze.html#evaluate_baseline_data

# https://doi.org/10.1111/eva.12787
# https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13209


library(assignPOP)

dat <- read.Structure("./data/6479snps_run4_K4structpops_modified.str", 
                      ploidy=2)

assign.MC( dat, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="svm", dir="./analysis/MC_svm_fst_propIndiv/")

assign.MC( dat, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="lda", dir="./analysis/MC_lda_fst_propIndiv/")

assign.MC( dat, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="naiveBayes", dir="./analysis/MC_naiveBayes_fst_propIndiv/")

assign.MC( dat, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="tree", dir="./analysis/MC_tree_fst_propIndiv/")

assign.MC( dat, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="randomForest", dir="./analysis/MC_RF_fst_propIndiv/")


# then run with fixed number of indivs, to see if we have bias bc of group size
assign.MC( dat, train.inds=c(5,10,12), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="svm", dir="./analysis/MC_svm_fst_fixedIndiv/")

assign.MC( dat, train.inds=c(5,10,12), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="lda", dir="./analysis/MC_lda_fst_fixedIndiv/")

assign.MC( dat, train.inds=c(5,10,12), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="naiveBayes", dir="./analysis/MC_naiveBayes_fst_fixedIndiv/")

assign.MC( dat, train.inds=c(5,10,12), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="tree", dir="./analysis/MC_tree_fst_fixedIndiv/")

assign.MC( dat, train.inds=c(5,10,12), train.loci=c(0.1, 0.25, 0.5, 0.75, 1),
           loci.sample="fst", iterations=30, model="randomForest", dir="./analysis/MC_RF_fst_fixedIndiv/")

# calculate accuracy:
accuracyMC_svm_fst_propIndiv <- accuracy.MC(dir= "./analysis/MC_svm_fst_propIndiv/")
accuracyMC_lda_fst_propIndiv <- accuracy.MC(dir= "./analysis/MC_lda_fst_propIndiv/")
accuracyMC_naiveBayes_fst_propIndiv <- accuracy.MC(dir= "./analysis/MC_naiveBayes_fst_propIndiv/")
accuracyMC_tree_fst_propIndiv <- accuracy.MC(dir= "./analysis/MC_tree_fst_propIndiv/")
accuracyMC_RF_fst_propIndiv <- accuracy.MC(dir= "./analysis/MC_RF_fst_propIndiv/")

accuracyMC_svm_fst_fixedIndiv <- accuracy.MC(dir= "./analysis/MC_svm_fst_fixedIndiv/")
accuracyMC_lda_fst_fixedIndiv <- accuracy.MC(dir= "./analysis/MC_lda_fst_fixedIndiv/")
accuracyMC_naiveBayes_fst_fixedIndiv <- accuracy.MC(dir= "./analysis/MC_naiveBayes_fst_fixedIndiv/")
accuracyMC_tree_fst_fixedIndiv <- accuracy.MC(dir= "./analysis/MC_tree_fst_fixedIndiv/")
accuracyMC_RF_fst_fixedIndiv <- accuracy.MC(dir= "./analysis/MC_RF_fst_fixedIndiv/")


# combine results

library(dplyr)

accuracyMC_svm_fst_propIndiv <- accuracyMC_svm_fst_propIndiv %>%
  mutate(test_set = "svm_propIndiv", .before = 1)
accuracyMC_lda_fst_propIndiv <- accuracyMC_lda_fst_propIndiv %>%
  mutate(test_set = "lda_propIndiv", .before = 1)
accuracyMC_naiveBayes_fst_propIndiv <- accuracyMC_naiveBayes_fst_propIndiv %>%
  mutate(test_set = "naiveBayes_propIndiv", .before = 1)
accuracyMC_tree_fst_propIndiv <- accuracyMC_tree_fst_propIndiv %>%
  mutate(test_set = "tree_propIndiv", .before = 1)
accuracyMC_RF_fst_propIndiv <- accuracyMC_RF_fst_propIndiv %>%
  mutate(test_set = "RF_propIndiv", .before = 1)

accuracyMC_svm_fst_fixedIndiv <- accuracyMC_svm_fst_fixedIndiv %>%
  mutate(test_set = "svm_fixedIndiv", .before = 1)
accuracyMC_lda_fst_fixedIndiv <- accuracyMC_lda_fst_fixedIndiv %>%
  mutate(test_set = "lda_fixedIndiv", .before = 1)
accuracyMC_naiveBayes_fst_fixedIndiv <- accuracyMC_naiveBayes_fst_fixedIndiv %>%
  mutate(test_set = "naiveBayes_fixedIndiv", .before = 1)
accuracyMC_tree_fst_fixedIndiv <- accuracyMC_tree_fst_fixedIndiv %>%
  mutate(test_set = "tree_fixedIndiv", .before = 1)
accuracyMC_RF_fst_fixedIndiv <- accuracyMC_RF_fst_fixedIndiv %>%
  mutate(test_set = "RF_fixedIndiv", .before = 1)

# combine dfs
propIndiv <- rbind(accuracyMC_svm_fst_propIndiv,accuracyMC_lda_fst_propIndiv, accuracyMC_naiveBayes_fst_propIndiv,
      accuracyMC_tree_fst_propIndiv, accuracyMC_RF_fst_propIndiv)

fixedIndiv <- rbind(accuracyMC_svm_fst_fixedIndiv, accuracyMC_lda_fst_fixedIndiv, accuracyMC_naiveBayes_fst_fixedIndiv,
      accuracyMC_tree_fst_fixedIndiv, accuracyMC_RF_fst_fixedIndiv)



# plot results:

library(ggplot2)
library(tidyr)
library(dplyr)

propIndiv$train.inds <- as.factor(propIndiv$train.inds)
propIndiv$train.loci <- as.factor(propIndiv$train.loci)

df <- propIndiv
pop <- c("all", "Orange_84", "Red_21", "DarkGreen_17","LightGreen_14" )
colnames(df) <- gsub("^assign\\.rate\\.", "", names(df))
#reshape data frame
df_melted <- df %>%
  pivot_longer(
    cols = c(all, Orange_84, DarkGreen_17, LightGreen_14, Red_21),
    names_to = "Population",
    values_to = "Assignment_Rate"
  )

head(df_melted)
levels(dfre$variable) <- sub("all", "Overall", levels(dfre$variable)) #change "all" to "Overall" if exists

plt_prop <- ggplot(df_melted, aes(x=train.inds, y=Assignment_Rate, fill=train.loci))+
  geom_boxplot()+
  facet_grid(test_set  ~ Population )+
  xlab("Proportion of individuals used in training set") + ylab("Assignment accuracy") +
  scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) 

ggsave("figures/assignment_rates_proportionIndivs.png",
       h=6,w=10)





# then plot the fixed indivs approach

df <- fixedIndiv

df$train.inds <- as.factor(df$train.inds)
df$train.loci <- as.factor(df$train.loci)

pop <- c("all", "Orange_84", "Red_21", "DarkGreen_17","LightGreen_14" )
colnames(df) <- gsub("^assign\\.rate\\.", "", names(df))
#reshape data frame
df_melted <- df %>%
  pivot_longer(
    cols = c(all, Orange_84, DarkGreen_17, LightGreen_14, Red_21),
    names_to = "Population",
    values_to = "Assignment_Rate"
  )

head(df_melted)
levels(dfre$variable) <- sub("all", "Overall", levels(dfre$variable)) #change "all" to "Overall" if exists
df_melted$train.inds <- factor(df_melted$train.inds, 
                               levels = c("5", "10", "12"))


plt_fixed <- ggplot(df_melted, aes(x=train.inds, y=Assignment_Rate, fill=train.loci))+
  geom_boxplot()+
  facet_grid(test_set  ~ Population )+
  xlab("Number of individuals used in training set") + ylab("Assignment accuracy") +
  scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) 

ggsave("figures/assignment_rates_fixedIndivs.png", plt_fixed,
       h=6,w=10)










#### --------------------------------------------------------------------------
#### --------------------------------------------------------------------------
#### --------------------------------------------------------------------------
#### --------------------------------------------------------------------------
#### --------------------------------------------------------------------------
# rubias
# https://cran.r-project.org/web/packages/rubias/vignettes/rubias-overview.html
library(rubias)
library(adegenet)

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

# change to rubias format.
df <- cbind(data.frame(sample_type = rep("reference", nrow(dat@tab)),
           repunit = as.character(dat@pop),
           collection = as.character(dat@pop),
           indiv = dat@other$X
           ),
           as.data.frame(dat@tab))

# self assign leave one out
self_df <- self_assign(reference = df, gen_start_col = 5)
head(self_df, n = 100)

# summarize repunit results
sa_to_repu <- self_df %>%
  group_by(indiv, collection, repunit) %>%
  top_n(1, scaled_likelihood) # just the top assignment for each sample

# summary of assignments without a likelihood threshold
assign_no_thres <- sa_to_repu %>%
  group_by(repunit, inferred_repunit) %>%
  tally()


# summarize assignments with 50% likelihood threshold:
# 50% likelihood threshold
thres50 <- self_df %>%
  group_by(indiv, collection, repunit) %>%
  filter(scaled_likelihood > 0.5) %>%
  group_by(repunit, inferred_repunit) %>%
  tally() %>%
  rename(threshold_50 = n)


# 90% likelihood threshold
thres90 <- self_df %>%
  group_by(indiv, collection, repunit) %>%
  filter(scaled_likelihood > 0.9) %>%
  group_by(repunit, inferred_repunit) %>%
  tally() %>%
  rename(threshold_90 = n)

# combine:
assign_no_thres %>%
  left_join(., thres50) %>%
  left_join(., thres90)



# summary of assignments with a likelihood threshold of 0.9
sa_to_repu %>%
  filter(scaled_likelihood > 0.9) %>%
  group_by(repunit, inferred_repunit) %>%
  tally() %>%
  ungroup() %>%
  group_by(repunit) %>%
  mutate(total = sum(n)) %>%
  mutate(correct = ifelse(repunit == inferred_repunit, n/total, 0)) %>%
  ungroup() %>%
  filter(repunit == inferred_repunit) 

# very bad assignment rates... 

sa_to_repu %>%
  ggplot(aes(x = z_score)) +
  geom_density() +
  facet_wrap(.~collection)

# make number of samples even:
set.seed(765) # need to set a seed to make this reproducible!

rubias_genos_14 <- df %>%
  group_by(collection) %>%
  sample_n(14, replace = FALSE) %>% 
  ungroup() # lowest pop has 14, so use this

# perform self-assignment on reduced number of colony samples
assign14 <- self_assign(reference = rubias_genos_14, gen_start_col = 5)



# summarize repunit results
top_assign14 <- assign14 %>%
  group_by(indiv, collection, repunit) %>%
  top_n(1, scaled_likelihood) # just the top assignment for each sample

# summary of assignments without a likelihood threshold
assign14_no_thres <- top_assign14 %>%
  group_by(repunit, inferred_repunit) %>%
  tally()

# summarize with 50% likelihood:
# 50% likelihood threshold
thres50_14samples <- top_assign14 %>%
  group_by(indiv, collection, repunit) %>%
  filter(scaled_likelihood > 0.5) %>%
  group_by(repunit, inferred_repunit) %>%
  tally() %>%
  rename(threshold_50 = n)

# 90% likelihood threshold
thres90_14samples <- top_assign14 %>%
  group_by(indiv, collection, repunit) %>%
  filter(scaled_likelihood > 0.9) %>%
  group_by(repunit, inferred_repunit) %>%
  tally() %>%
  rename(threshold_90 = n)

# combine
assign14_no_thres %>%
  left_join(., thres50_14samples) %>%
  left_join(., thres90_14samples)


# summary of assignments without a likelihood threshold
top_assign14 %>%
  #filter(scaled_likelihood > 0.9) %>%
  group_by(repunit, inferred_repunit) %>%
  tally() %>%
  ungroup() %>%
  group_by(repunit) %>%
  mutate(total = sum(n)) %>%
  mutate(correct = ifelse(repunit == inferred_repunit, n/total, 0)) %>%
  ungroup() %>%
  filter(repunit == inferred_repunit) #%>%










sum(self_df$repunit == self_df$inferred_repunit)/nrow(self_df)

sa_to_repu <- self_df %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood))

head(sa_to_repu)

inferred <- sa_to_repu[which(sa_to_repu$repu_scaled_like > 0.7),]

inferred[!inferred$repunit == inferred$inferred_repunit,]
sum(inferred$repunit == inferred$inferred_repunit)/136






# Simulated mixtures using a leave-one-out type of approach

df_sims <- assess_reference_loo(reference = df, 
                                  gen_start_col = 5, 
                                  reps = 5, 
                                  mixsize = 100)


tmp <- chin_sims_repu_top6 %>%
  mutate(repunit = ifelse(repunit %in% arep$repunit, repunit, "OTHER")) %>%
  group_by(iter, repunit) %>%
  summarise(true_repprop = sum(true_pi), 
            reprop_posterior_mean = sum(post_mean_pi),
            repu_n = sum(n)) %>%
  mutate(repu_n_prop = repu_n / sum(repu_n))




df_sims

dim(dat@tab)

length(dat$LocusName)
ncol(dat$DataMatrix)

dat <- read.Structure("./data/6479snps_run4_K4structpops_modified.str", 
                      ploidy=2)
