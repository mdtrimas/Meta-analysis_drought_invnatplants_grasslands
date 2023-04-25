################################################################################
########## Meta-analysis - analyzing extracted data ############

#for Mac
setwd("~/Box Sync/metaanalysis_droughtinvasion")
search()


#loading packages
install.packages("vegan")
install.packages("tidyverse")
install.packages("codyn")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("plotly")
install.packages("nmle")
install.packages("lme4")
install.packages("olsrr")
install.packages("car")
install.packages("patchwork")
install.packages("lmerTest")
install.packages("piecewiseSEM")
install.packages("multcomp")
install.packages("MuMIn")
install.packages("forcats")
install.packages("emmeans")

library(vegan)
library(tidyverse)
library(codyn)
library(ggplot2)
library(reshape2)
library(plotly)
library(nlme)
library(lme4)
library(olsrr)
library(car)
library(patchwork)
library(lmerTest)
library(piecewiseSEM)
library(multcomp)
library(MuMIn)
library(forcats)
library(emmeans)

#Set ggplot2 theme to black and white
theme_set(theme_bw())
#Update ggplot2 theme - make box around the x-axis title size 30, vertically justify x-axis title to 0.35, 
#Place a margin of 15 around the x-axis title.  
#Make the x-axis title size 30. For y-axis title, make the box size 30, put the writing at a 90 degree angle, and vertically justify the title to 0.5.  
#Add a margin of 15 and make the y-axis text size 25. Make the plot title size 30 and vertically justify it to 2.  Do not add any grid lines.  
#Do not add a legend title, and make the legend size 20
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=12)),
             axis.text.x=element_text(size=20), axis.title.y=element_text(size=20, angle=90, vjust=0.5,
                                                                          margin=margin(r=15)), axis.text.y=element_text(size=20), plot.title =
               element_text(size=20, vjust=2), panel.grid.major=element_blank(),
             panel.grid.minor=element_blank(),
             legend.text=element_text(size=15))


#set colorblind friendly color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000", "#CC79A7")



##################################################
####### Read in data ########
meta <- read.csv("data_extraction_meta.csv", header = TRUE, na.strings = "")


#number of unique studies = 66 studies
length(unique(meta$study_ID))





#################################################
######## Response ratios - Hedges' G ######
#calculate manually instead of using a package
meta_rr <- meta %>%
  mutate(control_n = as.numeric(control_n), drought_n = as.numeric(drought_n), driver_n = as.numeric(driver_n), drought_driver_n = as.numeric(drought_driver_n)) %>%
  #drought
  mutate(poolsd_dr = sqrt(((drought_n - 1) * drought_sd^2 + (control_n - 1) * control_sd^2) / (drought_n + control_n - 2))) %>% #pooled SD
  mutate(poolsd_dr = ifelse(poolsd_dr == 0, 0.00001, poolsd_dr)) %>%
  mutate(hedges_g_rr_dr = (drought_avg - control_avg) / poolsd_dr) %>%  #calculate RR
  mutate(resp_dir_dr = ifelse(hedges_g_rr_dr < 1, -1, 
                      ifelse(hedges_g_rr_dr > 1, 1, hedges_g_rr_dr))) %>%   #direction of response
  mutate(abs_hedges_dr = abs(hedges_g_rr_dr)) %>%    #absolute value of hedges
  mutate(ln_hedges_dr = log(abs_hedges_dr + 1)) %>% #lnrr + 1
  mutate(lnrr_dr = (resp_dir_dr * ln_hedges_dr)) %>%  #multiply response by direction 
  #driver
  mutate(poolsd_driv = sqrt(((driver_n - 1) * driver_sd^2 + (control_n - 1) * control_sd^2) / (driver_n + control_n - 2))) %>% #pooled SD
  mutate(poolsd_driv = ifelse(poolsd_driv == 0, 0.00001, poolsd_driv)) %>%
  mutate(hedges_g_rr_driv = (driver_avg - control_avg) / poolsd_driv) %>%  #calculate RR
  mutate(resp_dir_driv = ifelse(hedges_g_rr_driv < 1, -1, 
                              ifelse(hedges_g_rr_driv > 1, 1, hedges_g_rr_driv))) %>%   #direction of response
  mutate(abs_hedges_driv = abs(hedges_g_rr_driv)) %>%    #absolute value of hedges
  mutate(ln_hedges_driv = log(abs_hedges_driv + 1)) %>% #lnrr + 1
  mutate(lnrr_driv = (resp_dir_driv * ln_hedges_driv)) %>%  #multiply response by direction 
  #drought driver
  mutate(poolsd_drdriv = sqrt(((drought_driver_n - 1) * drought_driver_sd^2 + (control_n - 1) * control_sd^2) / (drought_driver_n + control_n - 2))) %>% #pooled SD
  mutate(poolsd_drdriv = ifelse(poolsd_drdriv == 0, 0.00001, poolsd_drdriv)) %>%
  mutate(hedges_g_rr_drdriv = (drought_driver_avg - control_avg) / poolsd_drdriv) %>%  #calculate RR
  mutate(resp_dir_drdriv = ifelse(hedges_g_rr_drdriv < 1, -1, 
                              ifelse(hedges_g_rr_drdriv > 1, 1, hedges_g_rr_drdriv))) %>%   #direction of response
  mutate(abs_hedges_drdriv = abs(hedges_g_rr_drdriv)) %>%    #absolute value of hedges
  mutate(ln_hedges_drdriv = log(abs_hedges_drdriv + 1)) %>% #lnrr + 1
  mutate(lnrr_drdriv = (resp_dir_drdriv * ln_hedges_drdriv))  #multiply response by direction 


#What are all the potential responses?
#aridity, inv_nat, funct_grp, study_type, measure_of_drought, mag_drt_dec, drought_length, resp_type




#################################################
###### Statistics ##########
#Overall lmer model
drought_lmer <- lmerTest::lmer(lnrr_dr ~ inv_nat + (1|study_ID), data = meta_rr)
anova(drought_lmer, type = 3)  #p = 0.0002891
AIC(drought_lmer)  #3274.914

drought_lmer2 <- lmerTest::lmer(lnrr_dr ~ inv_nat + (1|study_ID:scientific_name), data = meta_rr)
anova(drought_lmer2, type = 3)  #p = 0.02753
AIC(drought_lmer2) #3279.967

drought_lm <- lm(lnrr_dr ~ inv_nat, data = meta_rr)
anova(drought_lm)  #p = 0.003335
AIC(drought_lm)  #3288.434   - this has lowest AIC (simpler model), so I think we can go with simpler models and use glms



#big model
drought_lmer3 <- lmerTest::lmer(lnrr_dr ~ aridity + inv_nat + funct_grp + study_type + measure_of_drought + mag_drt_dec + drought_length + resp_type + (1|study_ID:scientific_name), data = meta_rr)
anova(drought_lmer3, type = 3)  #inv_nat p = 0.006401, study_type p = 0.053227
AIC(drought_lmer3)  #3341.762

drought_lmer4 <- lmerTest::lmer(lnrr_dr ~ aridity + inv_nat + funct_grp + study_type + measure_of_drought + mag_drt_dec + drought_length + resp_type + (1|study_ID), data = meta_rr)
anova(drought_lmer4, type = 3)  #only inv_nat is sig, p = 0.0002481
AIC(drought_lmer4)  #3339.435

drought_lm5 <- lm(lnrr_dr ~ aridity + inv_nat + funct_grp + study_type + measure_of_drought + mag_drt_dec + drought_length + resp_type, data = meta_rr)
anova(drought_lm5)  #aridity p = 0.0393944, inv_nat p = 0.0005071, funct_grp p = 0.0375227, study_type p = 0.0365013, mag_drt_dec = 0.0515412
AIC(drought_lm5)  #3280.915   - so far this has lowest AIC (simpler model)


#this one!
drought_lmer_overall <- lmerTest::lmer(lnrr_dr ~ inv_nat + (1|study_ID), data = meta_rr)
anova(drought_lmer_overall, type = 3)   #inv_nat p = 0.0002891


#aridity
drought_lmer_arid <- lmerTest::lmer(lnrr_dr ~ inv_nat * aridity + (1|study_ID), data = meta_rr)
anova(drought_lmer_arid, type = 3)   #inv_nat*arid p = 0.003773
AIC(drought_lmer_arid)  #3282.79   #maybe go with this bc i think this is the most technically correct, but can drop scientific name bc AIC is always higher
tuk_drought_arid <- emmeans(drought_lmer_arid, ~ inv_nat * aridity)
contrast(tuk_drought_arid, "consec", simple = "each", combine = TRUE, adjust = "mvt")  #within >0.65, n-i p = 0.0432; within 0.2-0.5, n-i p = 0.0015


#funct group
drought_lmer_funct <- lmerTest::lmer(lnrr_dr ~ inv_nat * funct_grp + (1|study_ID), data = meta_rr)
anova(drought_lmer_funct, type = 3)   #inv_nat*funct p = 0.08547

tuk_drought_funct <- emmeans(drought_lmer_funct, ~ inv_nat * funct_grp)
contrast(tuk_drought_funct, "consec", simple = "each", combine = TRUE, adjust = "mvt")  #within shrub, n-i p = 0.0196


#study type
drought_lmer_study <- lmerTest::lmer(lnrr_dr ~ inv_nat * study_type + (1|study_ID), data = meta_rr)
anova(drought_lmer_study, type = 3)   #inv_nat p = 0.06077

tuk_drought_study <- emmeans(drought_lmer_study, ~ inv_nat * study_type)
contrast(tuk_drought_study, "consec", simple = "each", combine = TRUE, adjust = "mvt")  #within greenshouse, n-i p = 0.0006


#measure of drought
drought_lmer_meas <- lmerTest::lmer(lnrr_dr ~ inv_nat * measure_of_drought + (1|study_ID), data = meta_rr)
anova(drought_lmer_meas, type = 3)   #inv_nat p = 0.0001835, inv_nat*measure p = 0.0909494

tuk_drought_meas <- emmeans(drought_lmer_meas, ~ inv_nat * measure_of_drought)
contrast(tuk_drought_meas, "consec", simple = "each", combine = TRUE, adjust = "mvt")  #within soil moist, n-i p = 0.00006


#magnitude of drought decrease
meta_rr_drdec <- meta_rr %>%
  drop_na(mag_drt_dec)
drought_lmer_dec <- lmerTest::lmer(lnrr_dr ~ inv_nat * mag_drt_dec + (1|study_ID), data = meta_rr_drdec)
anova(drought_lmer_dec, type = 3)   #inv_nat p = 0.01337

tuk_drought_dec <- emmeans(drought_lmer_dec, ~ inv_nat * mag_drt_dec)
contrast(tuk_drought_dec, "consec", simple = "each", combine = TRUE, adjust = "mvt")  #within 51-75, n-i p = 0.0142


#drought length
drought_lmer_len <- lmerTest::lmer(lnrr_dr ~ inv_nat * drought_length + (1|study_ID), data = meta_rr)
anova(drought_lmer_len, type = 3)   #inv_nat p = 0.02342; drought length p = 0.07979

tuk_drought_len <- emmeans(drought_lmer_len, ~ inv_nat * drought_length)
contrast(tuk_drought_len, "consec", simple = "each", combine = TRUE, adjust = "mvt")  #within 2-6 months, n-i p = 0.0049


#response type
drought_lmer_resp <- lmerTest::lmer(lnrr_dr ~ inv_nat * resp_type + (1|study_ID), data = meta_rr)
anova(drought_lmer_resp, type = 3)   #inv_nat p = 0.01118

tuk_drought_resp <- emmeans(drought_lmer_resp, ~ inv_nat * resp_type)
contrast(tuk_drought_resp, "consec", simple = "each", combine = TRUE, adjust = "mvt")  #within growth, n-i p = 0.0218

summary(glht(drought_lmer_resp, linfct = mcp(inv_nat = "Tukey"), test = adjusted(type = "BH"))) #all levels are very significantly different from one another

#response type 2
#combine abundance and production - use this one
meta_rr2 <- meta_rr %>%
  mutate(resp_type2 = ifelse(resp_type == "abundance", "production", resp_type))
drought_lmer_resp2 <- lmerTest::lmer(lnrr_dr ~ inv_nat * resp_type2 + (1|study_ID), data = meta_rr2)
anova(drought_lmer_resp2, type = 3)   #inv_nat p = 0.000232

#############################







############################################
######### Figures #########

cbPalette <- c("#0072B2", "#CC79A7", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#000000")

#overall inv vs native
meta_rr_inv_nat <- meta_rr %>%
  group_by(inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup()

inv_nat_fig <- ggplot(data = meta_rr_inv_nat, aes(x = avg, y = inv_nat, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.1, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "", y = "Plant Status", color = "Plant Status")


#response type
meta_rr_resp2 <- meta_rr2 %>%
  group_by(resp_type2, inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup() 

resp_fig2 <- ggplot(data = meta_rr_resp2, aes(x = avg, y = resp_type2, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "Effect Size", y = "Response Variable", color = "Plant Status")

inv_nat_fig / resp_fig2

#resp 800 x 800


#study type
meta_rr_study <- meta_rr %>%
  group_by(study_type, inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup()

study_fig <- ggplot(data = meta_rr_study, aes(x = avg, y = study_type, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "Effect Size", y = "Study Type", color = "Plant Status")



#aridity
meta_rr_arid <- meta_rr %>%
  group_by(aridity, inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup()

meta_rr_arid$aridity <- factor(meta_rr_arid$aridity, levels = c("0.03-0.2 arid", "0.2-0.5 semi-arid", "0.5-0.65 dry sub-humid", ">0.65 humid"))

arid_fig <- ggplot(data = meta_rr_arid, aes(x = avg, y = aridity, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1.2, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "", y = "Climate Class", color = "Plant Status")


#magnitude of drought decrease
meta_rr_dec <- meta_rr %>%
  drop_na(mag_drt_dec) %>%
  group_by(mag_drt_dec, inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup() 

dec_fig <- ggplot(data = meta_rr_dec, aes(x = avg, y = mag_drt_dec, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1.2, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "", y = "Percent Decrease in Water (%)", color = "Plant Status")


#drought length
meta_rr_length <- meta_rr %>%
  group_by(drought_length, inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup() 

meta_rr_length$drought_length <- factor(meta_rr_length$drought_length, levels = c("<2 months", "2-6 months", "6 months - 1.5 years", ">1.5 years"))

length_fig <- ggplot(data = meta_rr_length, aes(x = avg, y = drought_length, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1.2, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "Effect Size", y = "Duration of Drought", color = "Plant Status")

arid_fig + dec_fig + length_fig + plot_layout(ncol = 1)

#dr_fig 1200 x 1500

#measure of drought
meta_rr_meas <- meta_rr %>%
  group_by(measure_of_drought, inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup()

meas_fig <- ggplot(data = meta_rr_meas, aes(x = avg, y = measure_of_drought, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "Effect Size", y = "Drought Implementation", color = "Plant Status")



#functional group
meta_rr_func <- meta_rr %>%
  group_by(funct_grp, inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup()

func_fig <- ggplot(data = meta_rr_func, aes(x = avg, y = funct_grp, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-1, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "Effect Size", y = "Functional Group", color = "Plant Status")






#overall test for drought x driver
drought_lmer_overall_drdriv <- lmerTest::lmer(lnrr_drdriv ~ inv_nat + (1|study_ID), data = meta_rr)
anova(drought_lmer_overall_drdriv, type = 3)   #p = 0.5994


#overall test for drought (just those in the drought x driver group)
meta_rr_dr <- meta_rr %>%
  drop_na(driver_n)

drought_lmer_overall2 <- lmerTest::lmer(lnrr_dr ~ inv_nat + (1|study_ID), data = meta_rr_dr)
anova(drought_lmer_overall2, type = 3)   #inv_nat p = 0.001109



#overall inv vs native driver and droughtxdriver
meta_rr_inv_nat_dr <- meta_rr_dr %>%
  group_by(inv_nat) %>%
  summarise(std = sd(lnrr_dr, na.rm = TRUE), avg = mean(lnrr_dr, na.rm = TRUE), n = length(lnrr_dr)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup()

meta_rr_inv_nat_drdriv <- meta_rr %>%
  drop_na(lnrr_drdriv) %>%
  group_by(inv_nat) %>%
  summarise(std = sd(lnrr_drdriv, na.rm = TRUE), avg = mean(lnrr_drdriv, na.rm = TRUE), n = length(lnrr_drdriv)) %>%
  mutate(se = std/sqrt(n)) %>%
  mutate(marg = qt(0.975, df = n - 1) * avg/sqrt(n)) %>%
  mutate(lower = avg - marg) %>%
  mutate(upper = avg + marg) %>%
  ungroup()


inv_nat_fig_dr2 <- ggplot(data = meta_rr_inv_nat_dr, aes(x = avg, y = inv_nat, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-0.8, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "", y = "Plant Status", color = "Plant Status")


inv_nat_fig_drdriv <- ggplot(data = meta_rr_inv_nat_drdriv, aes(x = avg, y = inv_nat, color = inv_nat)) + 
  geom_vline(aes(xintercept = 0), size = 0.25, linetype = "dashed") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), linewidth = 0.8, height = 0.2, color = "gray50")+
  geom_point(size = 4) + 
  xlim(-0.8, 0.2) +
  scale_color_manual(values = cbPalette) + 
  labs(x = "Effect Size", y = "Plant Status", color = "Plant Status")


inv_nat_fig_dr2 + inv_nat_fig_drdriv + plot_layout(ncol = 1)
#dr_driv 800x800











###############################################
######### Map ########
# lots of libraries that mostly get the base map shapes and colors
install.packages("wesanderson")
devtools::install_github("ropenscilabs/rnaturalearth")
install.packages("rnaturalearthdata")
install.packages("sp")
install.packages("rgeos")
install.packages("maps")
install.packages("reshape")
install.packages("mapproj")

library(wesanderson)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(rgeos)
library(maps)
library(reshape)
library(mapproj)



#dataframe of only lat/long
meta2 <- meta %>%
  dplyr::select(c(study_ID, latitude, longitude))

meta3 <- meta2 %>%
  group_by(study_ID, latitude, longitude) %>%
  summarise(num = n()) %>%
  ungroup()



##this map
map <- ggplot(aes(size = num), data = meta3) +
  borders("world", colour = "gray40") +
  geom_point(data = meta3, mapping = aes(x = longitude, y = latitude), shape = 21, fill = "darkorchid4", alpha = 0.5) +  
  scale_size_continuous(range = c(4, 18)) +
  theme_bw() +
  theme(text = element_text(size = 24, colour = "black"), axis.text.x = element_text(size = 24, colour = "black"), axis.text.y = element_text(size = 24, colour = "black")) + #formatting the text
  ylab(expression("Latitude "*degree*"")) + #labels for the map x and y axes
  xlab(expression("Longitude "*degree*"")) +
  xlim(-170, 200) +
  ylim(-60, 80) + 
  labs(size = "Effect Size (#)") + #legend label
  theme(legend.position = c(0.15, 0.2), panel.grid.major = element_blank(), panel.grid.minor = element_blank())  #legend position

#with limits basically just cuts off Antarctica

#export at 1500 x 1000








