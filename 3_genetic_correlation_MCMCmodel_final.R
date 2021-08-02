# genetic correlation between traits  

# setwd("Z:/Anoob/MCMCglmm")  

## notes  
# gw - growth  
# bb - budbreak  
# bs - budset  
# B - lmer family blups (means)  
# P - plasticity value for the trait  

# notes from: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2017q1/025532.html  
# multiresponse mcmc models, page 92: https://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf  

# packages
require(MCMCglmm)
require(lmerTest)
require(dplyr)
require(tidyr)

# data 
budset_19 <- read.csv("./trait_data/BudSet_2019.csv")
budset_20 <- read.csv("./trait_data/BudSet_2020.csv")
budbreak_20 <- read.csv("./trait_data/BudBreak_2020_cGDD.csv")
heightGrowth_19 <- read.csv("./trait_data/Growth_2019.csv")
heightGrowth_20 <- read.csv("./trait_data/Growth_2020.csv")

## family blups
# blups for budset 2019
budset_19 <- budset_19 %>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)
budset_19 <- lmer(BudSet ~ Garden + (1|mBed) + (1|Family), budset_19)
budset_19 <- ranef(budset_19)
budset_19 <- budset_19$Family
budset_19 <- cbind(Family=rownames(budset_19),budset_19)
rownames(budset_19) <- NULL
names(budset_19)[2] <- "BudSet"
budset_19$Population <- gsub("\\_.*", "", budset_19$Family)

# blups for budset 2020
budset_20 <- budset_20 %>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)
budset_20 <- lmer(BudSet ~ Garden + (1|mBed) + (1|Family), budset_20)
budset_20 <- ranef(budset_20)
budset_20 <- budset_20$Family
budset_20 <- cbind(Family=rownames(budset_20),budset_20)
rownames(budset_20) <- NULL
names(budset_20)[2] <- "BudSet"
budset_20$Population <- gsub("\\_.*", "", budset_20$Family)

# blups for budbreak 2020
budbreak_20 <- budbreak_20 %>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)
budbreak_20 <- lmer(cGDD ~ Garden + (1|mBed) + (1|Family), budbreak_20)
budbreak_20 <- ranef(budbreak_20)
budbreak_20 <- budbreak_20$Family
budbreak_20 <- cbind(Family=rownames(budbreak_20),budbreak_20)
rownames(budbreak_20) <- NULL
names(budbreak_20)[2] <- "cGDD"
budbreak_20$Population <- gsub("\\_.*", "", budbreak_20$Family)

# blups for height growth 2019
heightGrowth_19 <- heightGrowth_19 %>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)
heightGrowth_19 <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_19)
heightGrowth_19 <- ranef(heightGrowth_19)
heightGrowth_19 <- heightGrowth_19$Family
heightGrowth_19 <- cbind(Family=rownames(heightGrowth_19),heightGrowth_19)
rownames(heightGrowth_19) <- NULL
names(heightGrowth_19)[2] <- "Growth"
heightGrowth_19$Population <- gsub("\\_.*", "", heightGrowth_19$Family)

# blups for height growth 2020
heightGrowth_20 <- heightGrowth_20 %>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)
heightGrowth_20 <- lmer(Growth ~ Garden + (1|mBed) + (1|Family), heightGrowth_20)
heightGrowth_20 <- ranef(heightGrowth_20)
heightGrowth_20 <- heightGrowth_20$Family
heightGrowth_20 <- cbind(Family=rownames(heightGrowth_20),heightGrowth_20)
rownames(heightGrowth_20) <- NULL
names(heightGrowth_20)[2] <- "Growth"
heightGrowth_20$Population <- gsub("\\_.*", "", heightGrowth_20$Family) 



# prior for genetic correlation for just one random effect 
prior_corr <- list(R = list(V = diag(2), nu = 0.002),
                   G = list(G1 = list(V = diag(2), nu = 2,
                                      alpha.mu = rep(0,2),
                                      alpha.V = diag(25^2,2,2))))

# Phenology vs growth----------------------------------------------------------------------------------#

######################################### Budset 2019-Growth 2019 ######################################

# against current year
bs_19_gw_19 <- merge(budset_19,heightGrowth_19[,c("Family","Growth")])
bs_19_gw_19 <- bs_19_gw_19[!is.na(bs_19_gw_19$BudSet),]
bs_19_gw_19 <- bs_19_gw_19[!is.na(bs_19_gw_19$Growth),]

# bs_19_gw_19 <- bs_19_gw_19%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bs_19_gw_19_fam_blup_pop_corr <- MCMCglmm(cbind(scale(BudSet),scale(Growth)) ~ trait - 1, 
                          random=~us(trait):Population,
                          rcov=~us(trait):units, # or idh
                          family=c("gaussian","gaussian"),
                          prior=prior_corr, 
                          # pedigree=Ped, 
                          data=bs_19_gw_19,
                          pr=TRUE, 
                          verbose=TRUE,
                          nitt=10000000, 
                          burnin=10000, 
                          thin=1000)

saveRDS(bs_19_gw_19_fam_blup_pop_corr,"./genetic_correlation_outputs/bs_19_gw_19_fam_blup_pop_corr")





######################################### Budset 2020-Growth 2020 ######################################

# against current year
bs_20_gw_20 <- merge(budset_20,heightGrowth_20[,c("Family","Growth")])
bs_20_gw_20 <- bs_20_gw_20[!is.na(bs_20_gw_20$BudSet),]
bs_20_gw_20 <- bs_20_gw_20[!is.na(bs_20_gw_20$Growth),]

# bs_20_gw_20 <- bs_20_gw_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bs_20_gw_20_fam_blup_pop_corr <- MCMCglmm(cbind(scale(BudSet),scale(Growth)) ~ trait - 1, 
                          random=~us(trait):Population,
                          rcov=~us(trait):units, # or idh
                          family=c("gaussian","gaussian"),
                          prior=prior_corr, 
                          # pedigree=Ped, 
                          data=bs_20_gw_20,
                          pr=TRUE, 
                          verbose=TRUE,
                          nitt=10000000, 
                          burnin=10000, 
                          thin=1000)

saveRDS(bs_20_gw_20_fam_blup_pop_corr,"./genetic_correlation_outputs/bs_20_gw_20_fam_blup_pop_corr")


######################################## Budbreak 2020-Growth 2019 #####################################

# against previous year
bb_20_gw_19 <- merge(budbreak_20,heightGrowth_19[,c("Family","Growth")])
bb_20_gw_19 <- bb_20_gw_19[!is.na(bb_20_gw_19$cGDD),]
bb_20_gw_19 <- bb_20_gw_19[!is.na(bb_20_gw_19$Growth),]

# bb_20_gw_19 <- bb_20_gw_19%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_gw_19_fam_blup_pop_corr <- MCMCglmm(cbind(scale(cGDD),scale(Growth)) ~ trait - 1, 
                          random=~us(trait):Population,
                          rcov=~us(trait):units, # or idh
                          family=c("gaussian","gaussian"),
                          prior=prior_corr, 
                          # pedigree=Ped, 
                          data=bb_20_gw_19,
                          pr=TRUE, 
                          verbose=TRUE,
                          nitt=10000000, 
                          burnin=10000, 
                          thin=1000)

saveRDS(bb_20_gw_19_fam_blup_pop_corr,"./genetic_correlation_outputs/bb_20_gw_19_fam_blup_pop_corr")






######################################## Budbreak 2020-Growth 2020 #####################################

# against current year
bb_20_gw_20 <- merge(budbreak_20,heightGrowth_20[,c("Family","Growth")])
bb_20_gw_20 <- bb_20_gw_20[!is.na(bb_20_gw_20$cGDD),]
bb_20_gw_20 <- bb_20_gw_20[!is.na(bb_20_gw_20$Growth),]

# bb_20_gw_20 <- bb_20_gw_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_gw_20_fam_blup_pop_corr <- MCMCglmm(cbind(scale(cGDD),scale(Growth)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bb_20_gw_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=1000)

saveRDS(bb_20_gw_20_fam_blup_pop_corr,"./genetic_correlation_outputs/bb_20_gw_20_fam_blup_pop_corr")



# Phenology vs phenology-------------------------------------------------------------------------------#

######################################## Budbreak 2020-Budset 2019 #####################################

# against previous year
bb_20_bs_19 <- merge(budbreak_20,budset_19[,c("Family","BudSet")])
bb_20_bs_19 <- bb_20_bs_19[!is.na(bb_20_bs_19$BudSet),]
bb_20_bs_19 <- bb_20_bs_19[!is.na(bb_20_bs_19$cGDD),]

# bb_20_bs_19 <- bb_20_bs_19%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_bs_19_fam_blup_pop_corr <- MCMCglmm(cbind(scale(cGDD),scale(BudSet)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bb_20_bs_19,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=1000)

saveRDS(bb_20_bs_19_fam_blup_pop_corr,"./genetic_correlation_outputs/bb_20_bs_19_fam_blup_pop_corr")






######################################## Budbreak 2020-Budset 2020 #####################################

# against current year
bb_20_bs_20 <- merge(budbreak_20,budset_20[,c("Family","BudSet")])
bb_20_bs_20 <- bb_20_bs_20[!is.na(bb_20_bs_20$BudSet),]
bb_20_bs_20 <- bb_20_bs_20[!is.na(bb_20_bs_20$cGDD),]

# bb_20_bs_20 <- bb_20_bs_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bb_20_bs_20_fam_blup_pop_corr <- MCMCglmm(cbind(scale(cGDD),scale(BudSet)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bb_20_bs_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=1000)

saveRDS(bb_20_bs_20_fam_blup_pop_corr,"./genetic_correlation_outputs/bb_20_bs_20_fam_blup_pop_corr")

# within trait correlation-----------------------------------------------------------------------------#
######################################### Budset 2019-Budset 2020 ######################################
budset_2019 <- budset_19
colnames(budset_2019)[colnames(budset_2019)=="BudSet"] <- "BudSet2019"
budset_2020 <- budset_20
colnames(budset_2020)[colnames(budset_2020)=="BudSet"] <- "BudSet2020"

bs_19_bs_20 <- merge(budset_2019,budset_2020[,c("Family","BudSet2020")])
bs_19_bs_20 <- bs_19_bs_20[!is.na(bs_19_bs_20$BudSet2020),]
bs_19_bs_20 <- bs_19_bs_20[!is.na(bs_19_bs_20$BudSet2019),]

# bs_19_bs_20 <- bs_19_bs_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

bs_19_bs_20_fam_blup_pop_corr <- MCMCglmm(cbind(scale(BudSet2019),scale(BudSet2020)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=bs_19_bs_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=1000)

saveRDS(bs_19_bs_20_fam_blup_pop_corr,"./genetic_correlation_outputs/bs_19_bs_20_fam_blup_pop_corr")

######################################### Growth 2019-Growth 2020 ######################################
heightGrowth_2019 <- heightGrowth_19
colnames(heightGrowth_2019)[colnames(heightGrowth_2019)=="Growth"] <- "Growth2019"
heightGrowth_2020 <- heightGrowth_20
colnames(heightGrowth_2020)[colnames(heightGrowth_2020)=="Growth"] <- "Growth2020"

gw_19_gw_20 <- merge(heightGrowth_2019,heightGrowth_2020[,c("Family","Growth2020")])
gw_19_gw_20 <- gw_19_gw_20[!is.na(gw_19_gw_20$Growth2019),]
gw_19_gw_20 <- gw_19_gw_20[!is.na(gw_19_gw_20$Growth2020),]

# gw_19_gw_20 <- gw_19_gw_20%>% unite(mBed, Garden, Bed, sep = "_", remove = FALSE)

gw_19_gw_20_fam_blup_pop_corr <- MCMCglmm(cbind(scale(Growth2019),scale(Growth2020)) ~ trait - 1, 
                             random=~us(trait):Population,
                             rcov=~us(trait):units, # or idh
                             family=c("gaussian","gaussian"),
                             prior=prior_corr, 
                             # pedigree=Ped, 
                             data=gw_19_gw_20,
                             pr=TRUE, 
                             verbose=TRUE,
                             nitt=10000000, 
                             burnin=10000, 
                             thin=1000)

saveRDS(gw_19_gw_20_fam_blup_pop_corr,"./genetic_correlation_outputs/gw_19_gw_20_fam_blup_pop_corr")

