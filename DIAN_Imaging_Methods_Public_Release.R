#### DIAN Imaging Methods
#### Code that accompanies the Nature Neuroscience DIAN Imaging Methods Manuscript
#### Author: Nicole S. McKay
#### April 2023

#### See the associated manuscript: XXXX
#### Please read the associated 'read me' file
#### This script assumes all data files are stored in a subfolder called 'Data'

#### Load packages ####
# assumes you already have the pacman package, if you do not: https://github.com/trinker/pacman 
#pacman::p_load(data.table, gt, gtsummary, effectsize, psych, tidyverse)
pacman::p_load(rstatix, broom, gt, gtsummary, psych, tidyverse, effectsize, ggalluvial, ggpubr)

#### loading each DIAN Data Release 15 file ####
# only selecting variables of interest, merging all data into one sheet
# creates a group variable to represent non-carriers, asymptomatic mutation carriers and symptomatic mutation carriers
# filters data to baseline visit only 
dian <- read.csv ("Data/DEMOGRAPHICS_DCF.csv")%>%
  select(newid15, visit, visitage, EDUC, SEX, HANDED, RACE, HISPANIC, RACEX, RACESEC, RACESECX, RACETER, RACETERX, PRIMLANG)

dian <- read.csv ("Data/GENETIC_DCF_delta_is_d.csv")%>%
  select(newid15, visit, master_famid, apoe, Mutation, MUTATIONTYPE, fam_mutation_chg)%>% 
  mutate(codon200 = case_when(Mutation == 1 & MUTATIONTYPE == 1 & parse_number(fam_mutation_chg) < 200 ~ "pre", Mutation == 1 & MUTATIONTYPE == 1 & fam_mutation_chg == "deletion exon 9" ~ "NA", Mutation == 1 & MUTATIONTYPE == 1 & parse_number(fam_mutation_chg) > 200 ~ "post"))%>%
  merge(.,dian, by = c("newid15", "visit"))

dian <- read.csv ("Data/psychometric_DCF.csv")%>%
  select(newid15, visit, MMSE,WAIS,MEMUNITS,ANIMALS, BOSTON)%>%
  merge(.,dian, by = c("newid15", "visit"))

dian <- read.csv ("Data/CLINICAL_DCF.csv")%>%
  select(newid15,visit, DIAN_EYO,cdrglob)%>%
  merge(.,dian, by = c("newid15", "visit"))

dian <- read.csv ("Data/IMAGING_dian_pj2_ADADCortSig.csv")%>%
  select(newid15, visit, CortSig_Thickness, MR_TOTV_HIPPOCAMPUS, MR_TOTV_INTRACRANIAL)%>%
  merge(.,dian, by = c("newid15", "visit"))

dian <- read.csv ("Data/FDG_DCF.csv")%>%
  select(newid15, visit, FDG_fSUVR_rsf_TOT_CTX_ISTHMUSCNG,FDG_fSUVR_rsf_TOT_CTX_INFERPRTL)%>%
  merge(.,dian, by = c("newid15", "visit"))

dian <- read.csv ("Data/pib_DCF.csv")%>%
  select(newid15, visit, PIB_fSUVR_rsf_TOT_CORTMEAN)%>%
  merge(.,dian, by = c("newid15", "visit"))%>%
  mutate(Group = case_when (Mutation == 1 & cdrglob == 0 ~ "MCA", Mutation == 1 & cdrglob > 0 ~ "MCS", Mutation == 0 ~ "NC"))%>%
  filter(visit == "v00")

#### Data Transforms #### 
individuals <- unique(dian$newid15) # 583 individuals who had imaging data collected

# remove individuals missing imaging data (ie their data failed QC), or missing genetic data 
dian <- dian%>% 
  drop_na(MR_TOTV_HIPPOCAMPUS, apoe)

# store information about flemmish/dutch mutations 
mutations <- dian%>% filter(fam_mutation_chg == "Glu693Gln" | fam_mutation_chg == "Ala692Gly")%>% mutate(count = "1")%>% select(fam_mutation_chg, count)

# run a regression to estimate influence of headsize specifically on hippocampal volume
hipp = lm(dian$MR_TOTV_HIPPOCAMPUS ~ dian$MR_TOTV_INTRACRANIAL); Bw = as.numeric(hipp$coefficients[2]); meICV = mean(dian$MR_TOTV_INTRACRANIAL) 

# apply the transform to hippocampal volumes, create FDG "summary" metric, and rename all imaging variables for easier future use 
dian <- dian %>%
  mutate(hippvol = MR_TOTV_HIPPOCAMPUS - (Bw * (MR_TOTV_INTRACRANIAL-meICV)), fdgsum = (FDG_fSUVR_rsf_TOT_CTX_ISTHMUSCNG + FDG_fSUVR_rsf_TOT_CTX_INFERPRTL)/2)%>%
  mutate(cortsig = CortSig_Thickness, pibsum = PIB_fSUVR_rsf_TOT_CORTMEAN)

# create a data frame of unimpaired controls with EYO between -10 and 0
controls <- dian%>%
  filter(Group == "NC" & cdrglob == 0.0 & DIAN_EYO <= 0 & DIAN_EYO >= -10)

# calculate z-scores for cognition and imaging variables of interest, relative to controls
mAnimals <- mean(controls$ANIMALS, na.rm = T); sdAnimals <- sd(controls$ANIMALS, na.rm = T); dian$zANIMALS <- (dian$ANIMALS - mAnimals)/sdAnimals
mMMSE <- mean(controls$MMSE, na.rm = T); sdMMSE <- sd(controls$MMSE, na.rm = T); dian$zMMSE <- (dian$MMSE - mMMSE)/sdMMSE


mWais <- mean(controls$WAIS, na.rm = T); sdWais <- sd(controls$WAIS, na.rm = T); dian$zWAIS <- (dian$WAIS - mWais)/sdWais
mMemunits <- mean(controls$MEMUNITS, na.rm = T); sdMemunits <- sd(controls$MEMUNITS, na.rm = T); dian$zMEMUNITS <- (dian$MEMUNITS - mMemunits)/sdMemunits
mBoston <- mean(controls$BOSTON, na.rm = T); sdBoston <- sd(controls$BOSTON, na.rm = T); dian$zBOSTON <- (dian$BOSTON - mBoston)/sdBoston
mPibsum <- mean(controls$pibsum, na.rm = T); sdPibsum <- sd(controls$pibsum, na.rm = T); dian$zpibsum <- (dian$pibsum - mPibsum)/sdPibsum
mFdgsum <- mean(controls$fdgsum, na.rm = T); sdFdgsum <- sd(controls$fdgsum, na.rm = T); dian$zfdgsum <- (dian$fdgsum - mFdgsum)/sdFdgsum
mCortsig <- mean(controls$cortsig, na.rm = T); sdCortsig <- sd(controls$cortsig, na.rm = T); dian$zcortsig <- (dian$cortsig - mCortsig)/sdCortsig
mHippvol <- mean(controls$hippvo, na.rm = T); sdHippvol <- sd(controls$hippvol, na.rm = T); dian$zhippvol <- (dian$hippvol - mHippvol)/sdHippvol
rm(mAnimals, sdAnimals, mWais, sdWais, mMemunits, sdMemunits, mBoston, sdBoston, mPibsum, sdPibsum, mFdgsum, sdFdgsum, mCortsig, sdCortsig, mHippvol, sdHippvol, Bw, meICV, hipp, mMMSE, sdMMSE)

# take the average of the z-scores for each person as the cognitive score
dian <- dian%>% 
  mutate(cognition = (zMEMUNITS + zWAIS + zANIMALS + zBOSTON)/4)%>% 
  mutate(Group = as.factor(Group), Group = fct_relevel(Group, c("NC", "MCA", "MCS")), SEX_bin = as.factor(SEX), 
         apoe = str_replace_all(apoe, c("34" = "Carrier", "33" = "Non-carrier", "44" = "Carrier", "23" = "Non-carrier", "24" = "Carrier", "22" = "Non-carrier")), 
         RACE_bin = case_when(RACE == "White" ~ "White", RACE != "White" ~ "Other"))

#### Demographics ####
individuals <- unique(dian$newid15)  # 534 individuals with imaging data 
families <- unique(dian$master_famid)  # 205 families 
unique_mutations <- unique(dian$fam_mutation_chg) # 108 different ADAD mutations represented
table(mutations) # 21 individuals with dutch mutation, 2 indiviudals with flemmish mutations

# Create a demographic table by participant group (substituting words in places of coded variables)
dian%>%
  select(Group, visitage, EDUC, SEX, HANDED, RACE, apoe, Mutation, MUTATIONTYPE, codon200, MMSE, DIAN_EYO, cdrglob, cognition)%>%
  mutate(SEX = str_replace_all(SEX, c("1" = "Male", "2" = "Female")), 
         Mutation = str_replace_all(Mutation, c("1" = "Carrier", "0" = "Non-carrier")), 
         MUTATIONTYPE = case_when(Mutation == "Carrier" & MUTATIONTYPE == "1" ~ "PSEN1",
                                  Mutation == "Carrier" & MUTATIONTYPE == "2" ~ "PSEN2", Mutation == "Carrier" & MUTATIONTYPE == "3" ~ "APP"), 
         cdrglob = case_when(cdrglob == 0.0 ~ "Unimpaired", cdrglob > 0.0 ~ "Impaired"), 
         Group = fct_relevel(Group, c("NC", "MCA", "MCS")))%>% 
  tbl_summary(by = Group, missing_text = "missing", type = all_continuous() ~ "continuous2", statistic = all_continuous() ~ c("{mean} ({sd})"))%>%
  bold_labels()%>% italicize_levels()%>% as_gt()
rm(individuals, families, mutations, unique_mutations)

#### Statistics ####
## Continuous variables - run anova with family as a covariate, calculate confidence intervals and effect sizes, perform follow up comparisons
# Age 
mod <- dian %>% anova_test(visitage ~ master_famid + Group); mod
confint((aov(visitage ~ Group + master_famid, data = dian)),level =0.95)
effectsize((aov(visitage ~ Group + master_famid, data = dian)), partial = TRUE)
mod2 <- dian %>% emmeans_test(visitage ~ Group, covariate = master_famid, p.adjust.method = 'bonferroni'); mod2
get_emmeans(mod2)

# Education 
mod <- dian %>% anova_test(EDUC ~ master_famid + Group); mod
confint((aov(EDUC ~ Group + master_famid, data = dian)),level =0.95)
effectsize((aov(EDUC ~ Group + master_famid, data = dian)), partial = TRUE)
mod2 <- dian %>% emmeans_test(EDUC ~ Group, covariate = master_famid, p.adjust.method = 'bonferroni'); mod2
get_emmeans(mod2)

# Categorical variables - run chi-square, calculate statistic and effect size 
# Sex
tbl = table(dian$SEX_bin, dian$Group)
mod <- chisq.test(tbl); tbl; mod
sqrt(mod$statistic / 534)

# Handedness 
tbl = table(dian$HANDED, dian$Group)
mod <- chisq.test(tbl); tbl; mod
sqrt(mod$statistic / 534)

# APOE e4
tbl = table(dian$apoe, dian$Group)
mod <- chisq.test(tbl); tbl; mod
sqrt(mod$statistic / 534)

# Race (binned by majority self-reported race)
tbl = table(dian$RACE_bin, dian$Group)
mod <- chisq.test(tbl); tbl; mod
sqrt(mod$statistic / 534)

# Since age and education differed across groups, we will use these as covariates for imaging and cognition analyses 

# Cognition
mod <- dian %>% anova_test(cognition ~ master_famid +visitage +SEX +Group); mod
confint((aov(cognition ~ Group + master_famid +visitage +SEX, data = dian)),level =0.95)
effectsize((aov(cognition ~ Group + master_famid +visitage +SEX, data = dian)), partial = TRUE)
mod2 <- dian %>% emmeans_test(cognition ~ Group, covariate = c(master_famid, visitage), p.adjust.method = 'bonferroni'); mod2
get_emmeans(mod2)

# Hippocampus
mod <- dian%>% drop_na(hippvol)%>% anova_test(hippvol ~ master_famid +visitage +SEX +Group); mod
confint((aov(hippvol ~ Group + master_famid +visitage +SEX, data = dian)),level =0.95)
effectsize((aov(hippvol ~ Group + master_famid +visitage +SEX, data = dian)), partial = TRUE)
mod2 <- dian %>% emmeans_test(hippvol ~ Group, covariate = c(master_famid, visitage), p.adjust.method = 'bonferroni'); mod2
get_emmeans(mod2)

# Cortical Signature
mod <- dian%>% drop_na(cortsig) %>% anova_test(cortsig ~ master_famid +visitage +SEX +Group); mod
confint((aov(cortsig ~ Group + master_famid +visitage +SEX, data = dian)),level =0.95)
effectsize((aov(cortsig ~ Group + master_famid +visitage +SEX, data = dian)), partial = TRUE)
mod2 <- dian %>% emmeans_test(cortsig ~ Group, covariate = c(master_famid, visitage), p.adjust.method = 'bonferroni'); mod2
get_emmeans(mod2)

# PIB Summary
mod <- dian %>% drop_na(pibsum) %>% anova_test(pibsum ~ master_famid +visitage +SEX +Group); mod
confint((aov(pibsum ~ Group + master_famid +visitage +SEX, data = dian)),level =0.95)
effectsize((aov(pibsum ~ Group + master_famid +visitage +SEX, data = dian)), partial = TRUE)
mod2 <- dian %>% emmeans_test(pibsum ~ Group, covariate = c(master_famid, visitage), p.adjust.method = 'bonferroni'); mod2
get_emmeans(mod2)

# FDG Summary
mod <- dian %>% drop_na(fdgsum)%>% anova_test(fdgsum ~ master_famid +visitage +SEX +Group); mod
confint((aov(fdgsum ~ Group + master_famid +visitage +SEX, data = dian)),level =0.95)
effectsize((aov(fdgsum ~ Group + master_famid +visitage +SEX, data = dian)), partial = TRUE)
mod2 <- dian %>% emmeans_test(fdgsum ~ Group, covariate = c(master_famid, visitage), p.adjust.method = 'bonferroni'); mod2
get_emmeans(mod2)

rm(tbl, mod2, mod)

#### Figure 1: Demographic Plot ####
dian%>% 
  select(Group, visitage, SEX)%>%
  mutate(visitage = case_when(visitage > 17 & visitage < 33 ~ '17-32', visitage > 32.9 & visitage < 43 ~ '33-42', visitage > 42.9 & visitage < 71 ~ '43-70'), 
         SEX = str_replace_all(SEX, c("1" = "Male", "2" = "Female")))%>% 
  ggplot(aes(axis1 = Group, axis2 = visitage, axis3 = SEX))+
  scale_x_discrete(limits = c("Mutation Status", "Age", "Sex")) +
  geom_flow(aes(fill = Group), stat = "alluvium", aes.bind = "alluvia", lode.guidance = "forward", decreasing = FALSE)+
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.5, decreasing = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum), family = "sans"), decreasing = FALSE, size = 2.4) +
  scale_fill_manual(values = c("white", "white", "white", "white", "white", "tan1", "firebrick", "steelblue"))+
  labs (y="Culmulative Frequency", x = "")+ 
  theme(panel.border = element_blank(), axis.line = element_line(),  
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size = 7), axis.title = element_text(size=7), legend.position = "none")
  
#### Figure 2: Imaging results ####
mu1 <- dian %>% drop_na(zpibsum)%>% group_by(Group)%>% mutate (grp.mean = mean(zpibsum), Factor = "zpibsum")%>% select(Group, grp.mean, Factor) %>% ungroup(); mu1 <- unique(mu1)
mu2 <- dian %>% drop_na(zfdgsum)%>% group_by(Group)%>% mutate (grp.mean = -(mean(zfdgsum)), Factor = "zfdgsum")%>% select(Group, grp.mean, Factor) %>% ungroup(); mu2 <- unique(mu2)
mu3 <- dian %>% drop_na(zhippvol)%>% group_by(Group)%>% mutate (grp.mean = -(mean(zhippvol)), Factor = "zhippvol")%>% select(Group, grp.mean, Factor)%>% ungroup(); mu3 <- unique(mu3)
mu4 <- dian %>% drop_na(zcortsig)%>% group_by(Group)%>% mutate (grp.mean = -(mean(zcortsig)), Factor = "zcortsig")%>% select(Group, grp.mean, Factor)%>% ungroup(); mu4 <- unique(mu4)
mu <- rbind (mu1, mu2, mu3, mu4) %>% mutate (Value = as.factor(Factor))

dian %>% 
  mutate(zhippvol = -zhippvol, zcortsig = -zcortsig, zfdgsum = -zfdgsum)%>% 
  gather (Value, Biomarker, c(zhippvol, zcortsig, zfdgsum, zpibsum)) %>%
  mutate(Value = factor(Value, levels = c("zpibsum", "zfdgsum", "zhippvol", "zcortsig")))%>%
  ggplot (aes(x = Biomarker, fill = Group))+
  geom_density(alpha = 0.4)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Group),linetype="dashed")+ 
  scale_color_manual(values = c("steelblue","tan1", "firebrick"), labels = c("Non-Carriers", "Mutation-Carriers: \nAsymptomatic", "Mutation-Carriers: \nSymptomatic"))+
  scale_fill_manual(values = c("steelblue","tan1", "firebrick"), labels = c("Non-Carriers", "Mutation-Carriers: \nAsymptomatic", "Mutation-Carriers: \nSymptomatic"))+
  labs (y="Density", x = "Accumulation of Pathology (z-score, scaled)")+
  theme(panel.border = element_blank(), axis.line = element_line(),  
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size = 7), axis.title = element_text(size=7), legend.position = "bottom", legend.text = element_text(size=5), 
        legend.title = element_blank())+
  facet_wrap(Value ~ ., ncol = 2, scales = "free", labeller = as_labeller(c(zhippvol = "Hippocampal Atrophy", zcortsig = "Cortical Thinning", zfdgsum = "Hypometabolism", zpibsum = "Amyloid Deposition")))+
  theme(strip.placement = "outside")+
  theme(strip.text = element_text(size = 7))+
  theme(legend.position = "bottom")


#### Figure 3: Can not be replicated as it is created from data that would unblind researchers to participants ####

#### Figure 4: EYO vs Age ####
EYO <- dian %>% 
  filter(Group != "NC")%>% 
  mutate(EYO = round(DIAN_EYO, digits = 0), zhippvol = -zhippvol, zcortsig = -zcortsig, zfdgsum = -zfdgsum, MMSE = - zMMSE, cognition = -cognition)%>%
  gather(Factor, Difference, c(MMSE, cognition, zcortsig, zhippvol, zpibsum, zfdgsum))%>% mutate(Factor = as.factor(Factor))%>%
  mutate(Group2 = case_when(
    Factor == "zpibsum" ~ "PET", Factor == "zfdgsum" ~ "PET", 
    Factor == "zhippvol" ~ "MRI", Factor == "zcortsig" ~ "MRI", 
    Factor == "MMSE" ~ "COG", Factor == "cognition" ~ "COG"), Group2 = fct_relevel(Group2, c("PET", "MRI", "COG")), Factor = fct_relevel(Factor, c("zpibsum", "zfdgsum", "zcortsig", "zhippvol", "MMSE", "cognition")))%>%
  drop_na(Difference)%>%
  ggplot(aes (x = EYO, y = Difference, color = Factor, linetype = Factor))+
  geom_smooth(method = "loess", se = FALSE)+ 
  geom_vline(xintercept=0.0,linetype = "dashed", color = "black")+
  xlim(-20,15)+
  scale_color_manual(values = c("#6A569F", "#6A569F", "#658A56", "#658A56",  "#E09F3E", "#E09F3E"), labels = c("Amyloid Deposition", "Hypometabolism", "Cortical Atrophy","Hippocampal Atrophy", "Clinical Symptoms", 
                                                                                                               "Cognitive Decline"))+ 
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed","solid", "dashed"), 
                        labels = c("Amyloid Deposition", "Hypometabolism", "Cortical Atrophy","Hippocampal Atrophy", "Clinical Symptoms", 
                                   "Cognitive Decline"))+
  labs (x ="Estimated Years to Symptom Onset (EYO)", y = "Degree of Pathology \n(Standardized Difference)", color = "", linetype = "") +
  theme(panel.border = element_blank(), axis.line = element_line(),  
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size = 7), axis.title = element_text(size=7), legend.position = "bottom", legend.text = element_text(size=5), 
        legend.title = element_blank())+
  facet_wrap(Group2 ~ ., ncol = 3, scales = "free", labeller = as_labeller(c(PET = "PET Measures", COG = "Cognitive Measures", MRI = "MRI Measures")))+
  theme(strip.placement = "outside")+
  theme(strip.text = element_text(size = 7))+
  theme(legend.position = "bottom")


Age <- dian %>% 
  filter(Group != "NC")%>% 
  mutate(EYO = round(DIAN_EYO, digits = 0), zhippvol = -zhippvol, zcortsig = -zcortsig, zfdgsum = -zfdgsum, MMSE = - zMMSE, cognition = -cognition)%>%
  gather(Factor, Difference, c(MMSE, cognition, zcortsig, zhippvol, zpibsum, zfdgsum))%>% mutate(Factor = as.factor(Factor))%>%
  mutate(Group2 = case_when(
    Factor == "zpibsum" ~ "PET", Factor == "zfdgsum" ~ "PET", 
    Factor == "zhippvol" ~ "MRI", Factor == "zcortsig" ~ "MRI", 
    Factor == "MMSE" ~ "COG", Factor == "cognition" ~ "COG"), Group2 = fct_relevel(Group2, c("PET", "MRI", "COG")), Factor = fct_relevel(Factor, c("zpibsum", "zfdgsum", "zcortsig", "zhippvol", "MMSE", "cognition")))%>%
  drop_na(Difference)%>%
  ggplot(aes (x = visitage, y = Difference, color = Factor, linetype = Factor))+
  geom_smooth(method = "loess", se = FALSE)+ 
  scale_color_manual(values = c("#6A569F", "#6A569F", "#658A56", "#658A56",  "#E09F3E", "#E09F3E"), labels = c("Amyloid Deposition", "Hypometabolism", "Cortical Atrophy","Hippocampal Atrophy", "Clinical Symptoms", 
                                                                                                               "Cognitive Decline"))+ 
  scale_linetype_manual(values=c("solid", "dashed","solid", "dashed","solid", "dashed"), 
                        labels = c("Amyloid Deposition", "Hypometabolism", "Cortical Atrophy","Hippocampal Atrophy", "Clinical Symptoms", 
                                   "Cognitive Decline"))+
  labs (x ="Age (years)", y = "Degree of Pathology \n(Standardized Difference)", color = "", linetype = "") +
  theme(panel.border = element_blank(), axis.line = element_line(),  
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size = 7), axis.title = element_text(size=7), legend.position = "bottom", legend.text = element_text(size=5), 
        legend.title = element_blank())+
  facet_wrap(Group2 ~ ., ncol = 3, scales = "free", labeller = as_labeller(c(PET = "PET Measures", COG = "Cognitive Measures", MRI = "MRI Measures")))+
  theme(strip.placement = "outside")+
  theme(strip.text = element_text(size = 7))+
  theme(legend.position = "bottom")

ggarrange(EYO, Age, nrow=2, labels = c("A", "B"), common.legend = TRUE, legend="bottom", align = "v")

rm(Age, EYO)

#### Figure 5: Can not be replicated as it is created from data that would unblind researchers to participants ####
#### Figure 6: Can not be replicated as it is created from data that would unblind researchers to participants ####
#### Figure ED1: Breakdown of Race ####
#### Figure ED2: Heterogeneity of ADAD mutations ####
mu1 <- dian %>% drop_na(zpibsum)%>% filter(Group == "MCS" & MUTATIONTYPE != "2")%>% group_by(MUTATIONTYPE)%>% mutate (grp.mean = mean(zpibsum), Factor = "zpibsum")%>% select(MUTATIONTYPE, grp.mean, Factor) %>% ungroup(); mu1 <- unique(mu1)
mu2 <- dian %>% drop_na(zfdgsum)%>% filter(Group == "MCS" & MUTATIONTYPE != "2")%>% group_by(MUTATIONTYPE)%>% mutate (grp.mean = -(mean(zfdgsum)), Factor = "zfdgsum")%>% select(MUTATIONTYPE, grp.mean, Factor) %>% ungroup(); mu2 <- unique(mu2)
mu3 <- dian %>% drop_na(zhippvol)%>% filter(Group == "MCS" & MUTATIONTYPE != "2")%>% group_by(MUTATIONTYPE)%>% mutate (grp.mean = -(mean(zhippvol)), Factor = "zhippvol")%>% select(MUTATIONTYPE, grp.mean, Factor)%>% ungroup(); mu3 <- unique(mu3)
mu4 <- dian %>% drop_na(zcortsig)%>% filter(Group == "MCS" & MUTATIONTYPE != "2")%>% group_by(MUTATIONTYPE)%>% mutate (grp.mean = -(mean(zcortsig)), Factor = "zcortsig")%>% select(MUTATIONTYPE, grp.mean, Factor)%>% ungroup(); mu4 <- unique(mu4)
mu <- rbind (mu1, mu2, mu3, mu4) %>% mutate (Value = as.factor(Factor), MUTATIONTYPE = str_replace_all(MUTATIONTYPE, c("1" = "PSEN1", "3" = "APP")))

dian %>% 
  filter(Group == "MCS" & MUTATIONTYPE != "2")%>%
  mutate(zhippvol = -zhippvol, zcortsig = -zcortsig, zfdgsum = -zfdgsum, MUTATIONTYPE = str_replace_all(MUTATIONTYPE, c("1" = "PSEN1", "3" = "APP")))%>%
  gather (Value, Biomarker, c(zhippvol, zcortsig, zfdgsum, zpibsum)) %>%
  mutate(Value = factor(Value, levels = c("zpibsum", "zfdgsum", "zhippvol", "zcortsig")))%>%
  ggplot (aes(x = Biomarker, fill = as.factor(MUTATIONTYPE)))+
  geom_density(alpha = 0.4)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=MUTATIONTYPE),linetype="dashed")+ 
  scale_color_manual(values = c("#673C4F", "#7698B3"), labels = c("PSEN1",  "APP"))+
  scale_fill_manual(values = c("#673C4F","#7698B3"), labels = c("PSEN1",  "APP"))+
  labs (y="Density", x = "Accumulation of Pathology (z-score, scaled)")+
  theme(panel.border = element_blank(), axis.line = element_line(),  
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size = 7), axis.title = element_text(size=7), legend.position = "bottom", legend.text = element_text(size=5), 
        legend.title = element_blank())+
  facet_wrap(Value ~ ., ncol = 2, scales = "free", labeller = as_labeller(c(zhippvol = "Hippocampal Atrophy", zcortsig = "Cortical Thinning", zfdgsum = "Hypometabolism", zpibsum = "Amyloid Deposition")))+
  theme(strip.placement = "outside")+
  theme(strip.text = element_text(size = 7))

rm (mu1, mu2, mu3, mu4, mu)
#### Figure ED3: Control descriptives ####
dian %>%
  filter(Group == "NC" & cdrglob == 0)%>% 
  mutate(visitage = case_when(visitage > 17 & visitage <= 30 ~ '18-30', visitage >= 31 & visitage <= 40 ~ '31-40', visitage >= 41 & visitage < 100 ~ '41-70'), 
         SEX = str_replace_all(SEX, c("1" = "Male", "2" = "Female")), apoe = as.factor(fct_relevel(apoe, c("Carrier", "Non-carrier"))))%>%
  ggplot(aes(axis1 = visitage, axis2 = SEX, axis3 = apoe))+
  scale_x_discrete(limits = c("Age", "Sex", "APOE Status")) +
  geom_flow(aes(fill = visitage), stat = "alluvium", aes.bind = "alluvia", lode.guidance = "forward")+
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum), family = "sans"), decreasing = FALSE, size = 2.4) +
  scale_fill_manual(values = c("#4D5359","#508484","#79C99E", "white", "white", "white", "white"))+
  labs (y="Culmulative Frequency", x = "")+
  theme(panel.border = element_blank(), axis.line = element_line(),  
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size = 7), axis.title = element_text(size=7), legend.position = "none")


#### Figure ED4: Can not be replicated as it is created from data that would unblind researchers to participants ####
#### Figure ED5: Can not be replicated as it is created from data that would unblind researchers to participants ####