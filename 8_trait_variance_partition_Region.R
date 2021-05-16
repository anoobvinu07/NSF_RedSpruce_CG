########################################################################################################
#------------------------------------trait variance partitioning---------------------------------------#
########################################################################################################

## notes
# gw - growth
# bb - budbreak
# bs - budset
# B - lmer family blups (means)
# P - plasticity value for the trait

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

# add mBed  
budset_19$mBed <- paste0(budset_19$Garden,"_",budset_19$Bed)
budset_20$mBed <- paste0(budset_20$Garden,"_",budset_20$Bed)
budbreak_20$mBed <- paste0(budbreak_20$Garden,"_",budbreak_20$Bed)
heightGrowth_19$mBed <- paste0(heightGrowth_19$Garden,"_",heightGrowth_19$Bed)
heightGrowth_20$mBed <- paste0(heightGrowth_20$Garden,"_",heightGrowth_20$Bed)

########################################################################################################
#-----------------------------------------------Priors-------------------------------------------------#
########################################################################################################
# set priors  
## testing prior [https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q4/021287.html]   
# current prior based on page 68 of the course notes for MCMCglmm  
# variance structure based on page 70, Compound Variance Structures  
# https://cran.r-project.org/web/packages/MCMCglmm/vignettes/CourseNotes.pdf  
Prior <- list(R=list(V = 1, nu = 0.002),           # R - prior on residual variance  
              G = list(                                  # G - list of priors for random effect variance  
                       G1 = list(V = diag(3) * 0.02, nu = 4),   
                       G2 = list(V = 1, nu = 0.002),
                       G3 = list(V = 1, nu = 0.002)))


## prior details for idh variance structures   
# An uninformative prior for the correlations is an improper prior with V=diag(dim(V))???0 and nu=dim(V)+1

# equivalent lmer model
# mod1 <- lmer(BudSet ~ (1|Garden) + (Garden|Family) + (1|mBed) , data=budset_19[budset_19$Region == "Core",])


############################################ Budset 2019 ###############################################

## Core
bs_19_mod_Core <-  MCMCglmm(BudSet ~ 1, 
                            random= ~ idh(1+Garden):Family + Garden + mBed,
                            family="gaussian",
                            data=budset_19[budset_19$Region=="Core",],
                            prior=Prior, 
                            pr=TRUE, # save random effects 
                            # burnin=100 , nitt=100000 , thin=100) # testing
                            burnin=10000 , nitt=10000000 , thin=1000)



saveRDS(bs_19_mod_Core,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bs_19_mod_Core")

## Margin
bs_19_mod_Margin <-  MCMCglmm(BudSet ~ 1, 
                              random= ~ idh(1+Garden):Family + Garden + mBed,
                              family="gaussian",
                              data=budset_19[budset_19$Region=="Margin",],
                              prior=Prior, 
                              pr=TRUE, 
                              burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bs_19_mod_Margin,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bs_19_mod_Margin")

## Edge
bs_19_mod_Edge <-  MCMCglmm(BudSet ~ 1, 
                            random= ~ idh(1+Garden):Family + Garden + mBed,
                            family="gaussian",
                            data=budset_19[budset_19$Region=="Edge",],
                            prior=Prior, 
                            pr=TRUE, 
                            burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bs_19_mod_Edge,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bs_19_mod_Edge")

############################################ Budset 2020 ###############################################

## Core
bs_20_mod_Core <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ idh(1+Garden):Family + Garden + mBed,
                          family="gaussian",
                          data=budset_20[budset_20$Region=="Core",],
                          prior=Prior, 
                          pr=TRUE, 
                          burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bs_20_mod_Core,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bs_20_mod_Core")

## Margin
bs_20_mod_Margin <-  MCMCglmm(BudSet ~ 1, 
                          random= ~ idh(1+Garden):Family + Garden + mBed,
                          family="gaussian",
                          data=budset_20[budset_20$Region=="Margin",],
                          prior=Prior, 
                          pr=TRUE, 
                          burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bs_20_mod_Margin,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bs_20_mod_Margin")

## Edge
bs_20_mod_Edge <- MCMCglmm(BudSet ~ 1, 
                         random= ~ idh(1+Garden):Family + Garden + mBed,
                         family="gaussian",
                         data=budset_20[budset_20$Region=="Edge",],
                         prior=Prior, 
                         pr=TRUE, 
                         burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bs_20_mod_Edge,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bs_20_mod_Edge")

########################################### Budbreak 2020 ##############################################

## Core
bb_20_mod_Core <-  MCMCglmm(cGDD ~ 1, 
                          random= ~ idh(1+Garden):Family + Garden + mBed,
                          family="gaussian",
                          data=budbreak_20[budbreak_20$Region=="Core",],
                          prior=Prior, 
                          pr=TRUE, 
                          burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bb_20_mod_Core,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bb_20_mod_Core")

## Margin
bb_20_mod_Margin <-  MCMCglmm(cGDD ~ 1, 
                          random= ~ idh(1+Garden):Family + Garden + mBed,
                          family="gaussian",
                          data=budbreak_20[budbreak_20$Region=="Margin",],
                          prior=Prior, 
                          pr=TRUE, 
                          burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bb_20_mod_Margin,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bb_20_mod_Margin")

## North Carolina
bb_20_mod_Edge <-  MCMCglmm(cGDD ~ 1, 
                            random= ~ idh(1+Garden):Family + Garden + mBed,
                            family="gaussian",
                            data=budbreak_20[budbreak_20$Region=="Edge",],
                            prior=Prior, 
                            pr=TRUE, 
                            burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(bb_20_mod_Edge,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/bb_20_mod_Edge")

############################################ Growth 2019 ###############################################

## Core
gw_19_mod_Core <-  MCMCglmm(Growth ~ 1, 
                            random= ~ idh(1+Garden):Family + Garden + mBed,
                            family="gaussian",
                            data=heightGrowth_19[heightGrowth_19$Region=="Core",],
                            prior=Prior, 
                            pr=TRUE, 
                            burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(gw_19_mod_Core,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/gw_19_mod_Core")

## Margin
gw_19_mod_Margin <-  MCMCglmm(Growth ~ 1, 
                              random= ~ idh(1+Garden):Family + Garden + mBed,
                              family="gaussian",
                              data=heightGrowth_19[heightGrowth_19$Region=="Margin",],
                              prior=Prior, 
                              pr=TRUE, 
                              burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(gw_19_mod_Margin,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/gw_19_mod_Margin")

## Edge
gw_19_mod_Edge <-  MCMCglmm(Growth ~ 1, 
                            random= ~ idh(1+Garden):Family + Garden + mBed,
                            family="gaussian",
                            data=heightGrowth_19[heightGrowth_19$Region=="Edge",],
                            prior=Prior, 
                            pr=TRUE, 
                            burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(gw_19_mod_Edge,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/gw_19_mod_Edge")

############################################ Growth 2020 ###############################################

## Core
gw_20_mod_Core <-  MCMCglmm(Growth ~ 1, 
                            random= ~ idh(1+Garden):Family + Garden + mBed,
                            family="gaussian",
                            data=heightGrowth_20[heightGrowth_20$Region=="Core",],
                            prior=Prior, 
                            pr=TRUE, 
                            burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(gw_20_mod_Core,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/gw_20_mod_Core")

## Maryland
gw_20_mod_Margin <-  MCMCglmm(Growth ~ 1, 
                              random= ~ idh(1+Garden):Family + Garden + mBed,
                              family="gaussian",
                              data=heightGrowth_20[heightGrowth_20$Region=="Margin",],
                              prior=Prior, 
                              pr=TRUE, 
                              burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(gw_20_mod_Margin,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/gw_20_mod_Margin")

## North Carolina
gw_20_mod_Edge <-  MCMCglmm(Growth ~ 1, 
                          random= ~ idh(1+Garden):Family + Garden + mBed,
                          family="gaussian",
                          data=heightGrowth_20[heightGrowth_20$Region=="Edge",],
                          prior=Prior, 
                          pr=TRUE, 
                          burnin=10000 , nitt=10000000 , thin=1000)

saveRDS(gw_20_mod_Edge,"/home/Anoob/mydata/Anoob/MCMCglmm/region_varPart_outputs/gw_20_mod_Edge")












