library(assignPOP)
library(ggplot2)
library(tidyr)
library(dplyr)


dat <- read.Structure("./data/6479snps_run4_K4structpops_modified.str", 
                      ploidy=2)

# run assignment tests. cycling over both proportion of training indivs and training loci
# run each model, then can check after to see which performs best.
# first running with proportion of total individuals in training set
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
#levels(df_melted$variable) <- sub("all", "Overall", levels(df_melted$variable)) #change "all" to "Overall" 

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
df_melted_fixed <- df %>%
  pivot_longer(
    cols = c(all, Orange_84, DarkGreen_17, LightGreen_14, Red_21),
    names_to = "Population",
    values_to = "Assignment_Rate"
  )

head(df_melted_fixed)
#levels(df_melted$variable) <- sub("all", "Overall", levels(df_melted$variable)) #change "all" to "Overall"
df_melted_fixed$train.inds <- factor(df_melted_fixed$train.inds, 
                               levels = c("5", "10", "12"))

plt_fixed <- ggplot(df_melted_fixed, aes(x=train.inds, y=Assignment_Rate, fill=train.loci))+
  geom_boxplot()+
  facet_grid(test_set  ~ Population )+
  xlab("Number of individuals used in training set") + ylab("Assignment accuracy") +
  scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) 

ggsave("figures/assignment_rates_fixedIndivs.png", plt_fixed,
       h=6,w=10)


#-------------------------------------------------------------------------------
# that's a lot to digest, so parse it down just to svm, which clearly performs best

svm1 <- df_melted[grep("svm",df_melted$test_set),]
svm2 <- df_melted_fixed[grep("svm",df_melted_fixed$test_set),]
svm <- rbind(svm1, svm2)

plt_svm <- ggplot(svm, aes(x=train.inds, y=Assignment_Rate, fill=train.loci))+
  geom_boxplot()+
  facet_grid(Population ~ test_set, scales="free" )+
  xlab("Number of individuals used in training set") + ylab("Assignment accuracy") +
  scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) 
plt_svm

ggsave("figures/assignment_rates_svm.png", plt_svm,
       h=8,w=6)




