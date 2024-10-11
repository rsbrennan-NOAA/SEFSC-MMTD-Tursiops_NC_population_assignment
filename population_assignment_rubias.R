
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

# summary of assignments without a likelihood threshold
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
# does this make sense?


# ---------------------
# cycle over the random subsampling of indivs to make sure these results are consistent. 
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

ggplot(out, aes(x=repunit, y=correct, fill=repunit))+
  geom_boxplot()+
  ylab("% correct assignments") +
  xlab("population")

ggsave("figures/rubias_subset14indivs_accuracy.png",
       h=4, w=5)

+
  #facet_grid(test_set  ~ Population )+
  xlab("Proportion of individuals used in training set") +
  scale_fill_discrete(name="Prop. of\ntrain loci",guide=guide_legend(reverse=TRUE))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor=element_blank()) 




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
