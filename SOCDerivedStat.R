#===============================================================================
#                               Tatra Mountains 
#                             Incubation experiment
#
# This script analyzes root amended soils:
#       - SOC-derived respiration (root-derived respiration is in different file)
#       - to run the scipt smoothly, all objects from RootFreeStat.R must be in
#         environment because the analysis recycle them - the point is to quantify
#         the priming effect and the root-free soils represent reference
#       - constant thermal regimes and freeze-thaw cycles are analyze separately 
################################################################################

#===============================================================================
#Libraries and ggplot theme are loaded from separate files
#===============================================================================
mytheme <- readRDS('ggtheme.rds')
source("packages.R")
#===============================================================================
# Data
#===============================================================================
TrLongerFull <- read.csv('TrLongerFull.csv')
################################################################################

#===============================================================================
# 1. Stable temperatures (-5, 0, 5) + SOC-derived respiration only [non-linear mixed effects]
# - nested models are designed and compared 
# - the best one is selected
#===============================================================================
# ======================= Data preparation 
TrSoil <- TrLongerFull |> filter(Origin != "Roots" & StabTemp == 'Yes')
TrSoil$Locality <- as.factor(TrSoil$Locality)
TrSoil$Origin <- as.factor(TrSoil$Origin)
# ---- Convert months of incubation to days of incubation to match SoilResp in ug/d/d ----
TrSoil$Time <- TrSoil$MonthsInc*30
nrow(TrSoil)
TrSoil_cc <- TrSoil[complete.cases(TrSoil[, c("SoilResp", "Temp", "Time", "Locality", "Origin")]), ]
nrow(TrSoil_cc)   # should be 288

# ======================= Non-linear model fitting (refiting root-free soil model) - nested models 
TrSoil_cc2 <- TrSoil_cc %>% left_join(toc_anchor_soil, by = "Locality")

# ---- Baseline: no root effect on either parameter, using N2's a-structure ----
Rc0 <- gnls(
  SoilResp ~ exp(logk) * plogis(logit_a) * TOC_anchor * exp(-exp(logk) * Time),
  data = TrSoil_cc2,
  params = list(logk ~ Temp, logit_a ~ Locality + Temp),
  start = c(logk = coef(N2)[c("logk.(Intercept)","logk.Temp")], logit_a = coef(N2)[grep("logit_a", names(coef(N2)))]),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)

# ---- Root addition affect the decomposition constant - true priming ---- 
Rc1 <- update(Rc0, params = list(logk ~ Temp + Origin, logit_a ~ Locality + Temp),
              start = c(coef(Rc0)[1:2], 0, coef(Rc0)[3:9]))
# ---- Root addition affect the decomposition constant differently at different temperatures - temperature dependent priming ---- 
Rc2 <- update(Rc1, params = list(logk ~ Temp * Origin, logit_a ~ Locality + Temp),
              start = c(coef(Rc1)[1:3], 0, coef(Rc1)[4:10]))
# ---- Root addition also affect fraction of available C - apparent priming ---- 
Rc3 <- update(Rc2, params = list(logk ~ Temp * Origin, logit_a ~ Locality + Temp + Origin),
              start = c(coef(Rc2), 0))

anova(Rc0, Rc1, Rc2, Rc3)
AIC(Rc0, Rc1, Rc2, Rc3)

#===============================================================================
# ------------------- Results of the best model Rc3---------------------------->
#===============================================================================
summary(Rc3)
R2M(Rc3)
anova(Rc3, type = 'marginal')

set.seed(123)
Rc3_boot <- boot_gnls(Rc3, nboot = 999)
boot_coefs_Rc3 <- as.data.frame(Rc3_boot$t)
names(boot_coefs_Rc3) <- names(Rc3_boot$t0)
boot_coefs_Rc3 <- boot_coefs_Rc3[complete.cases(boot_coefs_Rc3), ]
summary(boot_coefs_Rc3)

bT   <- boot_coefs_Rc3[["logk.Temp"]]
bR   <- boot_coefs_Rc3[["logk.OriginSoil with roots"]]
bTR  <- boot_coefs_Rc3[["logk.Temp:OriginSoil with roots"]]
ba   <- boot_coefs_Rc3[["logit_a.(Intercept)"]]
bAR  <- boot_coefs_Rc3[["logit_a.OriginSoil with roots"]]

# ---- effect of root addition on a ----
median((plogis(ba + bAR) - plogis(ba))*100)
quantile((plogis(ba + bAR) - plogis(ba))*100, probs = c(0.025, 0.975))

Q10_control <- exp(10*bT)
Q10_root    <- exp(10*(bT+bTR))

tibble(
  metric = c("Q10 control", "Q10 root-added", "logit_a shift from root (available fraction)"),
  median = c(median(Q10_control), median(Q10_root), median(bAR)),
  lwr = c(quantile(Q10_control,0.025), quantile(Q10_root,0.025), quantile(bAR,0.025)),
  upr = c(quantile(Q10_control,0.975), quantile(Q10_root,0.975), quantile(bAR,0.975))
)

# ---- visualization of priming effect ----
temp_grid <- seq(-5, 5, by = 0.5)

k_traj_Rc3 <- lapply(c("absent","present"), function(rt) {
  sapply(temp_grid, function(T) {
    k_draws <- k_Rc3fun(boot_coefs_Rc3, T, rt)
    c(median = median(k_draws), lwr = quantile(k_draws,0.025), upr = quantile(k_draws,0.975))
  }) %>% t() %>% as_tibble() %>% mutate(Temp = temp_grid, Root = rt)
}) %>% bind_rows()

names(k_traj_Rc3)[1:3] <- c("k_median","k_lwr","k_upr")

# ---- fold change ratio between -5 and +5 ----
a_neg5_roots <- plogis(boot_coefs_Rc3[["logit_a.(Intercept)"]] + boot_coefs_Rc3[["logit_a.LocalityRA"]] +
                         boot_coefs_Rc3[["logit_a.OriginSoil with roots"]] +
                         boot_coefs_Rc3[["logit_a.Temp"]]*-5)
a_5_roots    <- plogis(boot_coefs_Rc3[["logit_a.(Intercept)"]] + boot_coefs_Rc3[["logit_a.LocalityRA"]] +
                         boot_coefs_Rc3[["logit_a.OriginSoil with roots"]] +
                         boot_coefs_Rc3[["logit_a.Temp"]]*5)

a_ratio_roots <- a_5_roots / a_neg5_roots

tibble(
  metric = "Fold-change in available fraction, -5 to 5C",
  median = median(a_ratio_roots),
  lwr = quantile(a_ratio_roots, 0.025),
  upr = quantile(a_ratio_roots, 0.975)
)


#=======Soil differencess in alpha
for(i in unique(TrLongerFull$Locality)){
  
  if(i == 'FU'){
    alfa_i =  plogis(boot_coefs_Rc3[['logit_a.(Intercept)']] + boot_coefs_Rc3[['logit_a.OriginSoil with roots']])
  } else {
    alfa_i = plogis(boot_coefs_Rc3[['logit_a.(Intercept)']] + boot_coefs_Rc3[['logit_a.OriginSoil with roots']] +
                      boot_coefs_Rc3[[paste0('logit_a.Locality', i)]])
  }
  
  print(c(i, mean(alfa_i)*100))
}



# ---- Carbon loss estimate ----
SoilCLoss <- data.frame(Locality = rep(unique(TrLongerFull$Locality), each = 3),
                        Temp = rep(c(-5, 0, 5), times = 6),
                        Center = numeric(length = 6*3),
                        lw = numeric(length = 6*3),
                        up = numeric(length = 6*3))
SoilCLoss$Locality <- factor(SoilCLoss$Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA"))

for(i in unique(SoilCLoss$Locality)){
  for(j in c(-5, 0, 5)){
    closs <- a_Rc3fun(boot_coefs_Rc3, i, j, root = "present")*(1-exp(-k_Rc3fun(boot_coefs_Rc3, j, root = "present")*120))*100
    SoilCLoss[(SoilCLoss$Locality == i & SoilCLoss$Temp == j), 3:5] <- c(mean(closs),
                                                                         quantile(closs, probs = c(0.025, 0.975)))
  }
}
for(i in  c("SC", "VW", "NW", "FU", "JA", "RA")){
  print(SoilCLoss[SoilCLoss$Locality == i, ])
}

#===============================================================================
# 2.Freeze-thaw cycles (-5/5, -18/5) + SOC-derived respiration in root-amended soils only
# - the best model Rc3 is used to make predictions on respiration in 
#   alternating thermal regimes after 2 months of incubation
# - metanalytical approach is used to quantify difference between predictions and 
#   measurements
# - the trajectory of freeze-thaw effect is extrapolated to later incubations
#===============================================================================
# ---- The same pipeline as the one for root-free soils is implemented ----
# ---- Get the coefficients ----
k_Rc3fun <- function(draws, temp, root) {
  base <- draws[["logk.(Intercept)"]] + draws[["logk.Temp"]] * temp
  if (root == "present") {
    base <- base + draws[["logk.OriginSoil with roots"]] + draws[["logk.Temp:OriginSoil with roots"]] * temp
  }
  exp(base)
}

a_Rc3fun <- function(draws, locality, temp, root) {
  base <- draws[["logit_a.(Intercept)"]]
  off <- if (locality == localities[1]) 0 else draws[[paste0("logit_a.Locality", locality)]]
  base <- base + off + draws[["logit_a.Temp"]] * temp
  if (root == "present") {
    base <- base + draws[["logit_a.OriginSoil with roots"]]
  }
  plogis(base)
}

# ---- Predictions for month-2 ----
predict_month2_Rc3 <- function(locality, root) {
  toc_val <- toc_anchor_soil$TOC_anchor[toc_anchor_soil$Locality == locality]
  
  k_neg5 <- k_Rc3fun(boot_coefs_Rc3, -5, root)
  k_5    <- k_Rc3fun(boot_coefs_Rc3, 5, root)
  a_neg5 <- a_Rc3fun(boot_coefs_Rc3, locality, -5, root)
  a_5    <- a_Rc3fun(boot_coefs_Rc3, locality, 5, root)
  
  pred_stable_neg5 <- k_neg5 * a_neg5 * toc_val * exp(-k_neg5 * 60.88)
  pred_stable_5    <- k_5    * a_5    * toc_val * exp(-k_5 * 60.88)
  
  tibble(
    Series = c("stable_-5", "stable_5"),
    pred_mean = c(mean(pred_stable_neg5), mean(pred_stable_5)),
    pred_lwr  = c(quantile(pred_stable_neg5, 0.025), quantile(pred_stable_5, 0.025)),
    pred_upr  = c(quantile(pred_stable_neg5, 0.975), quantile(pred_stable_5, 0.975)),
    Locality = locality,
    Root = root
  )
}

pred_val <- bind_rows(
  lapply(localities, function(l) predict_month2_Rc3(l, "absent")),
  lapply(localities, function(l) predict_month2_Rc3(l, "present"))
)

print(pred_val, n=30)

# ---- Observations ----
TrSoil_cc2$Root <- ifelse(TrSoil_cc2$Origin == 'Control soil', 'absent', 'present')
unique(TrSoil_cc2$Root)

obs_month2_Rc3 <- TrSoil_cc2 %>%
  filter(Temp %in% c(-5, 5), round((Time - 60.88)/30.44) == 0) %>%
  mutate(Series = ifelse(Temp == -5, "stable_-5", "stable_5"),
         Root = as.character(Root)) %>%
  group_by(Series, Locality, Root) %>%
  summarize(obs_mean = mean(SoilResp, na.rm=TRUE), n = n(), .groups="drop")

print(obs_month2_Rc3, n=30)

# ---- validation check for stable temperatures ----
compare_val <- pred_val %>%
  left_join(obs_month2_Rc3, by = c("Series", "Locality", "Root")) %>%
  mutate(within_band = obs_mean >= pred_lwr & obs_mean <= pred_upr)

print(compare_val, n=30)

compare_val %>% group_by(Root) %>% summarize(n_within = sum(within_band), n_total = n())

# ---- quantifying freeze-thaw effect ----
# join against observed SoilResp from root-added FT jars, compute residuals, same as the control version
FT_rootadded_df <- TrLongerFull %>%
  filter(Origin == 'Soil with roots' & TempFactor %in% c("5", "-5", "-18/5", "-5/5")) %>%
  mutate(TempFactor = factor(TempFactor, levels = c("5", "-5", "-18/5", "-5/5")),
         Locality = factor(Locality))
FT_rootadded_df$Time <- FT_rootadded_df$MonthsInc*30

obs_FT_rootadded <- FT_rootadded_df %>%
  filter(TempFactor %in% c("5", "-5", "-5/5", "-18/5")) %>%
  mutate(Time_since_thaw = round((Time - 60.88) / 30.44)) %>%   # months since thaw, from days
  rename(Treatment = TempFactor)

predict_5_5_Rc3 <- function(locality, root) {
  toc_val <- toc_anchor_soil$TOC_anchor[toc_anchor_soil$Locality == locality]
  k_neg5 <- k_Rc3fun(boot_coefs_Rc3, -5, root)
  k_5    <- k_Rc3fun(boot_coefs_Rc3, 5, root)
  a_5    <- a_Rc3fun(boot_coefs_Rc3, locality, 5, root)
  
  pred <- k_5 * a_5 * toc_val * exp(-(k_neg5 * 60.88 + k_5 * 14))
  tibble(
    pred_mean = mean(pred),
    pred_lwr  = quantile(pred, 0.025),
    pred_upr  = quantile(pred, 0.975)
  )
}

pred_5_5_val <- bind_rows(
  lapply(localities, function(l) bind_cols(Locality=l, Root="absent",  predict_5_5_Rc3(l,"absent"))),
  lapply(localities, function(l) bind_cols(Locality=l, Root="present", predict_5_5_Rc3(l,"present")))
)

print(pred_5_5_val, n=30)

obs_5_5_Rc3 <- obs_FT_rootadded %>%
  filter(Treatment == "-5/5", Time_since_thaw == 0) %>%
  mutate(Root = "present") %>%
  bind_rows(
    obs_FT %>% filter(Treatment == "-5/5", Time_since_thaw == 0) %>% mutate(Root = "absent")
  ) %>%
  group_by(Locality, Root) %>%
  summarize(obs_mean = mean(SoilResp, na.rm=TRUE), n = n(), .groups="drop")

print(obs_5_5_Rc3, n=30)

# ---- rebuilding for -18/5 treatment (a(-18) is not necessary so we just assume close to 0 k at -18) ----
# ---- predictions ----
predict_18_5_Rc3 <- function(locality, root) {
  toc_val <- toc_anchor_soil$TOC_anchor[toc_anchor_soil$Locality == locality]
  k_5 <- k_Rc3fun(boot_coefs_Rc3, 5, root)
  a_5 <- a_Rc3fun(boot_coefs_Rc3, locality, 5, root)
  
  pred <- k_5 * a_5 * toc_val * exp(-k_5 * 14)   # no frozen-phase depletion term
  tibble(pred_mean = mean(pred), pred_lwr = quantile(pred,0.025), pred_upr = quantile(pred,0.975))
}

pred_18_5_val <- bind_rows(
  lapply(localities, function(l) bind_cols(Locality=l, Root="absent",  predict_18_5_Rc3(l,"absent"))),
  lapply(localities, function(l) bind_cols(Locality=l, Root="present", predict_18_5_Rc3(l,"present")))
)

print(pred_18_5_val, n=30)

# ---- observations ----
obs_18_5_Rc3 <- obs_FT_rootadded %>%
  filter(Treatment == "-18/5", Time_since_thaw == 0) %>%
  mutate(Root = "present") %>%
  bind_rows(
    obs_FT %>% filter(Treatment == "-18/5", Time_since_thaw == 0) %>% mutate(Root = "absent")
  ) %>%
  group_by(Locality, Root) %>%
  summarize(obs_mean = mean(SoilResp, na.rm=TRUE), .groups="drop")

# ---- alternative approach to account for variance in observations and predictions ----
obs_summary_root <- obs_FT_rootadded %>% filter(Treatment=="-5/5", Time_since_thaw==0) %>%
  group_by(Locality) %>%
  summarize(obs_mean=mean(SoilResp,na.rm=TRUE), obs_sd=sd(SoilResp,na.rm=TRUE), n=n(), .groups="drop")

meta_data_present <- pred_5_5_val %>% filter(Root=="present") %>%
  left_join(obs_summary_root, by="Locality") %>%
  mutate(obs_se = obs_sd/sqrt(n),
         pred_se = (pred_upr-pred_lwr)/(2*1.96),
         yi = log(obs_mean) - log(pred_mean),
         vi = (obs_se/obs_mean)^2 + (pred_se/pred_mean)^2)

meta_data_both <- bind_rows(
  meta_data %>% mutate(Root="absent"),        # from the control-soil rma above
  meta_data_present %>% mutate(Root="present")
) %>%
  mutate(Root = factor(Root, levels=c("absent","present")))

print(meta_data_both)

res_mv_simple <- rma.mv(yi, vi, mods = ~ Root, random = ~ 1 | Locality,
                        data = meta_data_both, method = "REML")

summary(res_mv_simple)

exp(c(control = res_mv_simple$b[1], root_shift = res_mv_simple$b[1] + res_mv_simple$b[2]))
exp(res_mv_simple$b[1] + res_mv_simple$b[2] + c(-1,1)*1.96*res_mv_simple$se[2])  # CI on the ratio itself

# -18/5 treatment 
obs_summary_root_18 <- obs_FT_rootadded %>% filter(Treatment=="-18/5", Time_since_thaw==0) %>%
  group_by(Locality) %>%
  summarize(obs_mean=mean(SoilResp,na.rm=TRUE), obs_sd=sd(SoilResp,na.rm=TRUE), n=n(), .groups="drop")

meta_data_present_18 <- pred_18_5_val %>% filter(Root=="present") %>%
  left_join(obs_summary_root_18, by="Locality") %>%
  mutate(obs_se = obs_sd/sqrt(n),
         pred_se = (pred_upr-pred_lwr)/(2*1.96),
         yi = log(obs_mean) - log(pred_mean),
         vi = (obs_se/obs_mean)^2 + (pred_se/pred_mean)^2)

meta_data_both_18 <- bind_rows(
  meta_data_18 %>% mutate(Root="absent"),        # from the control-soil rma above
  meta_data_present_18 %>% mutate(Root="present")
) %>%
  mutate(Root = factor(Root, levels=c("absent","present")))

print(meta_data_both_18)

res_mv_simple_18 <- rma.mv(yi, vi, mods = ~ Root, random = ~ 1 | Locality,
                           data = meta_data_both_18, method = "REML")
summary(res_mv_simple_18)

exp(c(control = res_mv_simple_18$b[1], root_shift = res_mv_simple_18$b[1] + res_mv_simple_18$b[2]))
exp(res_mv_simple_18$b[1] + res_mv_simple_18$b[2] + c(-1,1)*1.96*res_mv_simple_18$se[2])  # CI on the ratio itself

# ---- Direct comparison of measurements at month 2 ----
#1. -5/5
anova(aov(SoilResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 2 & TempFactor != '0' & TempFactor != '-18/5' & Origin == 'Soil with roots' & Roots == 'Yes')))
anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 2 & TempFactor != '0' & TempFactor != '-18/5'  & Origin == 'Soil with roots' & Roots == 'Yes')))

TukeyHSD(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
             subset = (MonthsInc == 2 & TempFactor != '0' & TempFactor != '-18/5'  & Origin == 'Soil with roots' & Roots == 'Yes')))

#1. -18/5
anova(aov(SoilResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 2 & TempFactor != '0' & TempFactor != '-5/5' & Origin == 'Soil with roots' & Roots == 'Yes')))
anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 2 & TempFactor != '0' & TempFactor != '-5/5'  & Origin == 'Soil with roots' & Roots == 'Yes')))

TukeyHSD(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
             subset = (MonthsInc == 2 & TempFactor != '0' & TempFactor != '-5/5'  & Origin == 'Soil with roots' & Roots == 'Yes')))


# ---- Check if the respiration excess continue or disapears over time ----
# ---- this analysis is performed on the relative changes in respiration rates in respect to reference at month 2 ----
rel_decline_root <- bind_rows(
  obs_FT_rootadded %>% filter(Treatment == "-5/5") %>% mutate(Series = "-5/5_root"),
  obs_FT_rootadded %>% filter(Treatment == "-18/5") %>% mutate(Series = "-18/5_root"),
  obs_FT_rootadded %>% filter(Treatment == "-5") %>% mutate(Series = "-5_root"),
  obs_FT_rootadded %>% filter(Treatment == "5") %>% mutate(Series = "5_root")
) %>%
  group_by(Series, Locality) %>%
  mutate(month2_val = mean(SoilResp[Time_since_thaw == 0], na.rm = TRUE)) %>%
  filter(!is.na(month2_val), month2_val > 0) %>%
  mutate(rel_change = SoilResp / month2_val) %>%
  ungroup()

rel_time_model_5 <- lm(rel_change ~ Treatment * Time_since_thaw, data = rel_decline_root %>% filter(Time_since_thaw > 0 & Roots == "Yes" & Treatment != "-18/5"))
anova(rel_time_model_5)
summary(rel_time_model_5)

rel_time_model_18 <- lm(rel_change ~ Treatment * Time_since_thaw, data = rel_decline_root %>% filter(Time_since_thaw > 0 & Roots == "Yes" & Treatment != "-5/5"))
anova(rel_time_model_18)
summary(rel_time_model_18)
#==============RA - month 2 -5/5 treatment is removed, because the outlier pushes everything in incorrect direction
rel_decline_root_clean <- rel_decline_root %>%
  group_by(Series, Locality) %>%
  filter(first(month2_val) > 2) %>%
  ungroup()

rel_decline_root_clean |> group_by(Treatment, Time_since_thaw) |> summarise(mean(rel_change, na.rm = T))

rel_time_model_5 <- lm(rel_change ~ Treatment * Time_since_thaw, data = rel_decline_root_clean %>% filter(Time_since_thaw > 0 & Roots == "Yes" & Treatment != "-18/5"))
anova(rel_time_model_5)
summary(rel_time_model_5)

rel_time_model_18 <- lm(rel_change ~ Treatment * Time_since_thaw, data = rel_decline_root_clean %>% filter(Time_since_thaw > 0 & Roots == "Yes" & Treatment != "-5/5"))
anova(rel_time_model_18)
summary(rel_time_model_18)


#==============time to excess respiration attenuation
set.seed(123)
n_boot <- 2000

# 1. Bootstrap draws of the month-0 excess ratio, from the rma() model's own sampling distribution
ratio0_draws_soc <- exp(rnorm(n_boot, mean = sum(coef(res_mv_simple)), sd = sqrt(sum(res_mv_simple$se^2))))

# 2. Bootstrap the relative-decline trajectory by resampling replicates within each Series/Time cell
rel_boot_soc <- replicate(n_boot, {
  boot_sample <- rel_decline_root_clean %>% filter(Series == "-5/5_root") %>%
    group_by(Locality, Time_since_thaw) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    summarize(mean_rel = mean(rel_change, na.rm=TRUE), .groups="drop")
  boot_sample$mean_rel
}) %>% t()   # n_boot x 3 matrix, columns = Time_since_thaw 0,1,2

colnames(rel_boot_soc) <- unlist(lapply(unique(TrLongerFull$Locality)[-c(4, 5)], paste0, c("-m0", "-m2", "-m3")))


# 3. Combine: ratio(m) = ratio(0) * rel_obs(m)/rel_obs(0), but rel_obs(0) = 1 by construction
ratio_traj_boot_soc <- ratio0_draws_soc * rel_boot_soc   # element-wise: each row is one bootstrap replicate's full trajectory

ratio_traj_summary_soc <- tibble(
  Legend = colnames(rel_boot_soc),
  median = apply(ratio_traj_boot_soc, 2, median, na.rm = T),
  lwr = apply(ratio_traj_boot_soc, 2, quantile, 0.025, na.rm = T),
  upr = apply(ratio_traj_boot_soc, 2, quantile, 0.975, na.rm = T)
)
print(ratio_traj_summary_soc)

# boot_coefs_N2 and ratio0_draws/rel_boot have different lengths -> resample to align
# (treats the two models' bootstraps as independent)
set.seed(123)

idx <- sample(nrow(boot_coefs_Rc3), n_boot, replace = TRUE)

kn_soc <- k_Rc3fun(boot_coefs_Rc3, -5, 'present')[idx]
kp_soc <- k_Rc3fun(boot_coefs_Rc3, 5, 'present')[idx]


scen_soc <- list(all_cold = kn_soc,
                 all_warm = kp_soc)

traj_soc <- lapply(scen_soc, function(e) ratio0_draws_soc * rel_boot_soc / cbind(1, exp(-e*30), exp(-e*30*2),
                                                                                 1, exp(-e*30), exp(-e*30*2),
                                                                                 1, exp(-e*30), exp(-e*30*2),
                                                                                 1, exp(-e*30), exp(-e*30*2)))

trajDFrame_soc <- data.frame(Scenario = rep(c("Cold", "Warm"), each = 2000*12),
                             Locality = rep(rep(unique(TrLongerFull$Locality)[-c(4, 5)], each = 3*2000), 2),
                             Month = rep(rep(c(0, 1, 2), each = 2000), times = 4*2),
                             Rm = as.numeric(c(traj_soc$all_cold, traj_soc$all_warm))) 

trajDFrame_soc |> group_by(Scenario, Month) |> filter(Rm < 1000) |> 
  summarise(Center = mean(Rm, na.rm = T),
            lb = quantile(Rm, probs = c(0.025), na.rm = T),
            ub = quantile(Rm, probs = c(0.975), na.rm = T))

#==============time to excess respiration attenuation
set.seed(123)
n_boot <- 2000

# 1. Bootstrap draws of the month-0 excess ratio, from the rma() model's own sampling distribution
ratio0_draws_soc18 <- exp(rnorm(n_boot, mean = sum(coef(res_mv_simple_18)), sd = sqrt(sum(res_mv_simple_18$se^2))))

# 2. Bootstrap the relative-decline trajectory by resampling replicates within each Series/Time cell
rel_boot_soc18 <- replicate(n_boot, {
  boot_sample <- rel_decline_root_clean %>% filter(Series == "-18/5_root") %>%
    group_by(Locality, Time_since_thaw) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    summarize(mean_rel = mean(rel_change, na.rm=TRUE), .groups="drop")
  boot_sample$mean_rel
}) %>% t()   # n_boot x 3 matrix, columns = Time_since_thaw 0,1,2

colnames(rel_boot_soc18) <- unlist(lapply(unique(TrLongerFull$Locality)[-c(5)], paste0, c("-m0", "-m2", "-m3")))


# 3. Combine: ratio(m) = ratio(0) * rel_obs(m)/rel_obs(0), but rel_obs(0) = 1 by construction
ratio_traj_boot_soc18 <- ratio0_draws_soc18 * rel_boot_soc18   # element-wise: each row is one bootstrap replicate's full trajectory

ratio_traj_summary_soc18 <- tibble(
  Legend = colnames(rel_boot_soc18),
  median = apply(ratio_traj_boot_soc18, 2, median, na.rm = T),
  lwr = apply(ratio_traj_boot_soc18, 2, quantile, 0.025, na.rm = T),
  upr = apply(ratio_traj_boot_soc18, 2, quantile, 0.975, na.rm = T)
)
print(ratio_traj_summary_soc18)

# boot_coefs_N2 and ratio0_draws/rel_boot have different lengths -> resample to align
scen_soc18 <- list(all_cold = rep(0, 2000),
                   all_warm = kp_soc)

traj_soc18 <- lapply(scen_soc18, function(e) ratio0_draws_soc18 * rel_boot_soc18 / cbind(1, exp(-e*30), exp(-e*30*2),
                                                                                         1, exp(-e*30), exp(-e*30*2),
                                                                                         1, exp(-e*30), exp(-e*30*2),
                                                                                         1, exp(-e*30), exp(-e*30*2),
                                                                                         1, exp(-e*30), exp(-e*30*2)))

trajDFrame_soc18 <- data.frame(Scenario = rep(c("Cold", "Warm"), each = 2000*15),
                               Locality = rep(rep(unique(TrLongerFull$Locality)[-c(5)], each = 3*2000), 2),
                               Month = rep(rep(c(0, 1, 2), each = 2000), times = 5*2),
                               Rm = as.numeric(c(traj_soc18$all_cold, traj_soc18$all_warm))) 

trajDFrame_soc18 |> group_by(Scenario, Month) |> filter(Rm < 1000) |> 
  summarise(Center = mean(Rm, na.rm = T),
            lb = quantile(Rm, probs = c(0.025), na.rm = T),
            ub = quantile(Rm, probs = c(0.975), na.rm = T))

# ---- Direct comparison of measurements at month 4 ----
#1. -5/5
anova(aov(SoilResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Origin == 'Soil with roots')))
anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Origin == 'Soil with roots')))
TukeyHSD(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
             subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Origin == 'Soil with roots')))
SOC5end <- as.data.frame(TukeyHSD(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
                                      subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Origin == 'Soil with roots')))$`TempFactor:Locality`)
SOC5end$comp <- rownames(SOC5end)

for(i in unique(TrLonger$Locality)){
  #Soil specific comparison
  comp5negative <- paste0('-5:', i, '--18/5:', i)
  comp5positive <- paste0('5:', i, '--18/5:', i)
  
  print(rbind(
    SOC5end[grepl(comp5negative, SOC5end$comp), c("diff", 'p adj', 'comp')],
    SOC5end[grepl(comp5positive, SOC5end$comp), c("diff", 'p adj', 'comp')]
  )[2:3, ])
  
}


#1. -18/5
anova(aov(SoilResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Origin == 'Soil with roots')))
anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Origin == 'Soil with roots')))
TukeyHSD(aov(SoilResp ~ TempFactor, data = TrLongerFull,
             subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Origin == 'Soil with roots')))



#===============================================================================
# Figures for root-amended soils - SOC-derived respiration
#===============================================================================
#=======================Visualizing the true priming effect (constant k)  
TruePriming <- data.frame(Legend = rep(c("Control", "Root-added"), each = nrow(boot_coefs_Rc3)),
                          k = c(exp(boot_coefs_Rc3[, 'logk.(Intercept)']),
                                exp(boot_coefs_Rc3[, 'logk.(Intercept)'] + 
                                      boot_coefs_Rc3[, 'logk.OriginSoil with roots']))) 

TruePriming_Fig <- grid.arrange(
  TruePriming |> group_by(Legend) |> summarise(y = median(k),
                                               lwr = quantile(k, 0.025),
                                               upr = quantile(k, 0.975)) |>
    ggplot(aes(Legend, y)) + geom_point(cex = 6,  pch = 21, aes(fill = Legend), show.legend = F) +
    geom_errorbar(aes(ymin = lwr, ymax = upr, color = Legend), width = 0.1, show.legend = F) +
    scale_fill_manual(values = c("grey30", "firebrick")) +
    scale_color_manual(values = c("grey30", "firebrick")) +
    ylim(0, 0.05) + scale_x_discrete(labels = c('Root-free', 'Root-amended')) +
    labs(x = " ", y = expression(paste("Decay rate constant (", day^{-1}, ")")), title = "A)") +
    mytheme + theme(legend.title = element_blank(), legend.position = c(0.25, 0.85)),
  ggplot(k_traj_Rc3, aes(Temp, k_median, color = Root, fill = Root)) +
    geom_ribbon(aes(ymin = k_lwr, ymax = k_upr), alpha = 0.25, color = NA) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c(absent = "grey30", present = "firebrick"),
                       labels = c("Root-free", "Root-amended")) +
    scale_fill_manual(values = c(absent = "grey30", present = "firebrick"), guide = "none") +
    labs(x = "Temperature (\u00b0C)", y = expression(paste("Decay rate constant (", day^{-1}, ")")),
         title = "B)") +
    mytheme + theme(legend.title = element_blank(), legend.position = c(0.25, 0.85),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank()) +
    ylim(0, 0.05), ncol = 2
)

ggsave("StatOut/TruePriming_Fig.pdf", TruePriming_Fig,
       width = 280, height = 120, units = "mm", device = cairo_pdf)

#==================Measured and predicted SOC-derived respiration rates in root-amended soils 
#Measured
RootAmendedStable <- TrSoil_cc2 |> filter(Origin == 'Soil with roots') |> group_by(Locality, Temp, Time) |> summarise(y = mean(SoilResp, na.rm = T),
                                                                                                                      ySD = sd(SoilResp, na.rm = T),
                                                                                                                      TOC_anchor = median(TOC_anchor, na.rm = T))
#Predicted
RootAmendedStable$yhat <- numeric(length = nrow(RootAmendedStable))
RootAmendedStable$lwr <- numeric(length = nrow(RootAmendedStable))
RootAmendedStable$upr <- numeric(length = nrow(RootAmendedStable))

names(boot_coefs_Rc3)


for(i in 1:nrow(RootAmendedStable)){
  #Model parameters for row i
  loc_i <- RootAmendedStable$Locality[i]
  k_i = exp(boot_coefs_Rc3[['logk.(Intercept)']] + boot_coefs_Rc3[['logk.OriginSoil with roots']] + 
              (boot_coefs_Rc3[['logk.Temp']] +
                 boot_coefs_Rc3[['logk.Temp:OriginSoil with roots']])*RootAmendedStable$Temp[i])
  if(loc_i == 'FU'){
    alfa_i =  plogis(boot_coefs_Rc3[['logit_a.(Intercept)']] + boot_coefs_Rc3[['logit_a.OriginSoil with roots']] +
                       boot_coefs_Rc3[['logit_a.Temp']]*RootAmendedStable$Temp[i])
  } else {
    alfa_i = plogis(boot_coefs_Rc3[['logit_a.(Intercept)']] + boot_coefs_Rc3[['logit_a.OriginSoil with roots']] +
                      boot_coefs_Rc3[[paste0('logit_a.Locality', loc_i)]] + 
                      boot_coefs_Rc3[['logit_a.Temp']]*RootAmendedStable$Temp[i])
  }
  
  #Prediction for row i
  RootAmendedStable[i, c('yhat', 'lwr', 'upr')] <- matrix(c(mean(RootAmendedStable$TOC_anchor[i]*
                                                                   k_i*alfa_i*exp(-k_i*RootAmendedStable$Time[i])), 
                                                            quantile(x = RootAmendedStable$TOC_anchor[i]*
                                                                       k_i*alfa_i*exp(-k_i*RootAmendedStable$Time[i]),
                                                                     probs = c(0.025, 0.975))), nrow = 1)
  
  
}

RootAmendedStable_SOCderivedFig <- RootAmendedStable %>%
  mutate(Locality = factor(Locality,
                           levels = c("SC", "VW", "NW", "FU", "JA", "RA"),
                           labels = c(expression(atop("2160 m", "scree")),
                                      expression(atop("2160 m", "meadow")),
                                      expression(atop("2130 m", "meadow")),
                                      expression(atop("1780 m", "meadow")),
                                      expression(atop("1470 m", "forest")),
                                      expression(atop("1330 m", "forest"))))) |> 
  ggplot(aes(Time/30, y)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 1), aes(fill = factor(Temp)),
           color = 'black') + facet_wrap(. ~ Locality,
                                         nrow = 2, labeller = label_parsed) +
  mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD, color = factor(Temp)), 
                          position = position_dodge(width = 1), show.legend = F) +
  scale_fill_manual(values = c('white', 'grey', 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  labs(color = 'Thermal regime', fill = "Thermal regime", x = 'Time (months)',
       y = expression(paste("SOC-derived respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
  new_scale_fill() +
  new_scale_color() +
  geom_point(cex = 4, pch = 21, aes(y = yhat, fill = factor(Temp)), position = position_dodge(width = 1),
             show.legend = F) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr, colour = factor(Temp)), 
                position = position_dodge(width = 1), width = 0.1, show.legend = F) +
  scale_fill_manual(values = c('red', 'red', 'red')) +
  scale_color_manual(values = c('red', 'red', 'red')) +
  theme(legend.position = c(0.08, 0.88))

ggsave("StatOut/RootAmendedStable_SOCderivedFig.pdf", RootAmendedStable_SOCderivedFig,
       width = 300, height = 180, units = "mm", device = cairo_pdf)

#==================Trajectory of freeze-thaw effect on SOC-derived respiration rates in root-amended soils 
trajDFrame_soc <- rbind(trajDFrame_soc, data.frame(Scenario = rep(c("Cold", "Warm"), each = 2*3),
                                                   Locality = rep(rep(unique(TrLongerFull$Locality)[c(4, 5)], each = 3), 2),
                                                   Month = rep(c(0, 1, 2), times = 4),
                                                   Rm = rep(NA, times = 12))) 
trajDFrame_soc$Series <- c('-5/5')
trajDFrame_soc18 <- rbind(trajDFrame_soc18, data.frame(Scenario = rep(c("Cold", "Warm"), each = 3),
                                                       Locality = rep(rep(unique(TrLongerFull$Locality)[c(5)], each = 3), 1),
                                                       Month = rep(c(0, 1, 2), times = 2),
                                                       Rm = rep(NA, times = 6))) 
trajDFrame_soc18$Series <- c('-18/5')

SOCDerivedFT_temporalFig <- rbind(trajDFrame_soc, trajDFrame_soc18) |> group_by(Series, Scenario, Locality, Month, .drop = F) |> 
  filter(Rm < 1000) |> 
  summarise(y = mean(Rm, na.rm = T),
            lb = quantile(Rm, probs = c(0.025), na.rm = T),
            ub = quantile(Rm, probs = c(0.975), na.rm = T)) |> 
  mutate(Time = as.factor(Month + 2)) |> 
  ggplot(aes(factor(Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA")), y)) + 
  geom_point(cex = 6, pch = 21, aes(fill = Time), position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymin = lb, ymax = ub, color = Time), width = 0.1, , position = position_dodge(width = 1)) +
  facet_wrap(factor(Series, levels = c("-5/5", "-18/5"),
                    labels = c("-5˚C ⟶ +5˚C", "-18˚C ⟶ +5˚C")) ~ Scenario, scales = "free_y") + mytheme +
  geom_hline(yintercept = 1, color = 'red', lwd = 1.1) +
  scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                              expression(atop("2160 m", "meadow")),
                              expression(atop("2130 m", "meadow")),
                              expression(atop("1780 m", "meadow")),
                              expression(atop("1470 m", "forest")),
                              expression(atop("1330 m", "forest")))) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 14)) +
  labs(y = c("Relative effect of alternating thermal regime (nondimensional)"),
       color = c("Time (months)"), fill = c("Time (months)"))

ggsave("StatOut/SOCDerivedFT_temporalFig.pdf", SOCDerivedFT_temporalFig,
       width = 360, height = 180, units = "mm", device = cairo_pdf)

#=======================Freeze-thaw cycles - composite figure
RootAmendedRef <- TrSoil_cc |>  filter(Origin == "Soil with roots") |> group_by(Locality, Temp, MonthsInc) |> 
  summarise(yref = mean(SoilResp, na.rm = T), yrefSD = sd(SoilResp, na.rm = T)) |> 
  filter(Temp != 0) |> filter(MonthsInc != 3) |> 
  mutate(TempF = factor(Temp, levels = c(-5, 5), labels = c("-5˚C", '+5˚C')))

pred_FT_all_2 <- rbind(pred_FT_all[pred_FT_all$Root == 'present', ], pred_FT_all[pred_FT_all$Root == 'present', ])
pred_FT_all_2 <- pred_FT_all_2[, -2]
pred_FT_all_2$MonthsInc <- rep(c(2, 4), each = 12) 
pred_FT_all_2[13:nrow(pred_FT_all_2), c(2:4)] <- NA

SOCDerivedFT_compositeFig <- TrLongerFull |> filter(Origin == 'Soil with roots' & StabTemp == 'No' & MonthsInc != 3) |> 
  group_by(Locality, TempFactor, MonthsInc) |> summarise(y = mean(SoilResp, na.rm = T),
                                                         ySD = sd(SoilResp, na.rm = T)) |> 
  inner_join(pred_FT_all_2, by = c('Locality', "TempFactor", "MonthsInc")) |> 
  full_join(RootAmendedRef) |> 
  mutate(Month = factor(MonthsInc, levels = c(2, 4), labels = c("Month 2", "Month 4"))) |> 
  ggplot(aes(factor(Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA")), y)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 1),
           color = 'black', fill = 'white') +
  mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD), 
                          position = position_dodge(width = 1), show.legend = F) +
  #Reference constant temperature
  geom_bar(aes(y = yref, fill = TempF), stat = 'identity', width = 0.5,
           position = position_dodge(width = 0.5), alpha = 0.5) +
  geom_errorbar(aes(ymin = yref, ymax = yref + yrefSD, fill = TempF), stat = 'identity', width = 0.01,
                position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c('grey70', 'grey30')) +
  #=================================================================
facet_grid(Month ~ factor(TempFactor, levels = c("-5/5", "-18/5"),
                          labels = c("-5˚C ⟶ +5˚C", "-18˚C ⟶ +5˚C"))) +
  labs(color = 'Temperature regime', fill = "Constant thermal regime",
       y = expression(paste("SOC-derived respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
  new_scale_fill() +
  new_scale_color() +
  geom_point(cex = 2, pch = 21, aes(y = pred_mean, fill = TempFactor), position = position_dodge(width = 1),
             show.legend = F) + 
  geom_errorbar(aes(ymin = pred_lwr, ymax = pred_upr, colour = TempFactor), 
                position = position_dodge(width = 1), width = 0.1, show.legend = F, linewidth = 0.4) +
  scale_fill_manual(values = c('red', 'red')) +
  scale_color_manual(values = c('red', 'red')) +
  theme(axis.title.x = element_blank(), legend.position = c(0.15, 0.9), axis.text.x = element_text(size = 16)) +
  scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                              expression(atop("2160 m", "meadow")),
                              expression(atop("2130 m", "meadow")),
                              expression(atop("1780 m", "meadow")),
                              expression(atop("1470 m", "forest")),
                              expression(atop("1330 m", "forest"))))

ggsave("StatOut/SOCDerivedFT_compositeFig.pdf", SOCDerivedFT_compositeFig,
       width = 320, height = 180, units = "mm", device = cairo_pdf)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF SCRIPT
