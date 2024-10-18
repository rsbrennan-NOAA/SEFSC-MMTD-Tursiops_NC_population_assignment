
# rubias
# https://cran.r-project.org/web/packages/rubias/vignettes/rubias-overview.html
library(rubias)
library(adegenet)
library(dplyr)
library(ggplot2)
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

# change to rubias format.
df <- cbind(data.frame(sample_type = rep("reference", nrow(dat@tab)),
                       repunit = as.character(dat@pop),
                       collection = as.character(dat@pop),
                       indiv = dat@other$X),
            as.data.frame(dat@tab)
            )

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

# summary of assignments with a likelihood threshold of 0.9
assign_with_thres <- sa_to_repu %>%
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
assign_with_thres

write.csv(file="analysis/rubias_assign_rates_allIndivs.csv",assign_with_thres, 
          row.names=F)

# this is a known bias (see https://doi.org/10.1111/2041-210X.14286)


p1 <- ggplot(assign_with_thres, aes(x=repunit, y=correct, fill=repunit))+
  geom_point(size=4, pch=21)+
  ylab("% correct assignments") +
  xlab("population") + 
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none") +
  ylim(0,1) +
  scale_fill_manual(values=c("darkgreen", "lawngreen", "orange3", "red3")) +
  ggtitle("all individuals")

ggsave("figures/rubias_allindivs_accuracy.png",p1,
       h=4, w=5)

#--------------------------------------------------------------------------------
# make number of samples even:
set.seed(765) # need to set a seed to make reproducible

rubias_genos_14 <- df %>%
  group_by(collection) %>%
  sample_n(14, replace = FALSE) %>% 
  ungroup() # lowest pop has 14, so use this

# perform self-assignment on reduced number of samples
assign14 <- self_assign(reference = rubias_genos_14, gen_start_col = 5)

# summarize repunit results
top_assign14 <- assign14 %>%
  group_by(indiv, collection, repunit) %>%
  top_n(1, scaled_likelihood) # just the top assignment for each sample

# summary of assignments without a likelihood threshold
assign14_no_thres <- top_assign14 %>%
  group_by(repunit, inferred_repunit) %>%
  tally()

# summary of assignments with a likelihood threshold
top_assign14 %>%
  filter(scaled_likelihood > 0.9) %>%
  group_by(repunit, inferred_repunit) %>%
  tally() %>%
  ungroup() %>%
  group_by(repunit) %>%
  mutate(total = 14) %>%
  #mutate(total = sum(n)) %>%
  mutate(correct = ifelse(repunit == inferred_repunit, n/14, 0)) %>%
  ungroup() %>%
  filter(repunit == inferred_repunit) 

# seems to vastly improve with even sample sizes.



# -------------------------------------------------------------------------------
# cycle over the random subsampling of indivs to make sure these results are consistent w/ diff indivs. 
nreps <- 100

out <- data.frame(matrix(ncol=5, nrow=0))
colnames(out) <- c("repunit", "inferred_repunit", "n", "total", "correct")

for(i in 1:nreps){

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
  
  # count % accurate
  tmpout <- top_assign14 %>%
    filter(scaled_likelihood > 0.9) %>%
    group_by(repunit, inferred_repunit) %>%
    tally() %>%
    ungroup() %>%
    group_by(repunit) %>%
    mutate(total = 14) %>%
    #mutate(total = sum(n)) %>%
    mutate(correct = ifelse(repunit == inferred_repunit, n/14, 0)) %>%
    ungroup() %>%
    filter(repunit == inferred_repunit) 
  
  out <- rbind(out, tmpout)
  
  
  }


# plot results

p2 <- ggplot(out, aes(x=repunit, y=correct, fill=repunit))+
  geom_boxplot()+
  ylab("% correct assignments") +
  xlab("population") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position="none") +
  ylim(0,1) +
  scale_fill_manual(values=c("darkgreen", "lawngreen", "orange3", "red3")) +
  ggtitle("subsampled to 14 indivs")

ggsave("figures/rubias_subset14indivs_accuracy.png",p2,
       h=4, w=5)

ggsave("figures/rubias_combined_accuracy.png",ggarrange(p1, p2), h=3.5, w=7)
