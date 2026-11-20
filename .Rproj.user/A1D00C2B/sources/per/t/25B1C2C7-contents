#===============================================================================
#                               Tatra Mountains 
#                             Incubation experiment
#
# This script analyzes root-free soils:
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
# 1. Stable temperatures (-5, 0, 5) + soil respiration only [non-linear mixed effects]
# - nested models are designed and compared 
# - the best one is selected
#===============================================================================
# ======================= Data preparation 
Ctr <- TrLongerFull |> filter(Origin == "Control soil" & StabTemp == 'Yes')
Ctr$Locality <- as.factor(Ctr$Locality)

# ---- Per-soil per-temperature lm fits to guess starting values ----
start_vals_lin <- Ctr %>%
  group_by(Locality, Temp) %>%
  group_modify(~{
    d <- .x
    fit <- lm(log(SoilResp) ~ Time, data = d)
    k_est  <- -coef(fit)["Time"]
    C0_est <- exp(coef(fit)["(Intercept)"]) / k_est
    tibble(k = unname(k_est), C0 = unname(C0_est))
  })

print(start_vals_lin)

# ---- Data frame needs to be properly arranged so that the joint nonlinear model can be fit ----
Ctr <- Ctr |> arrange(Origin, MonthsInc, Locality, Temp)
# ---- Convert soil organic C from mg/g to ug/g ----
Ctr$soil_Cug <- Ctr$soil_C*1000
# ---- Convert months of incubation to days of incubation to match SoilResp in ug/d/d ----
Ctr$Time <- Ctr$MonthsInc*30
# ---- Creating unique identifier for accuracy control ----
Ctr$triplicate <- rep(seq(1,3), times = 6*3*3) 
Ctr$jar_id <- paste0(Ctr$Locality, '-', Ctr$Time, '-', Ctr$Temp, '-', Ctr$triplicate)

# ======================= Non-linear model fitting - first guesses
# ---- starting estimates of nonlinear mixed effect model build on first order decay rate ----
start_vals_lin_clean <- start_vals_lin %>%
  mutate(
    k_clean = ifelse(k <= 0, NA, k),
    C0_clean = ifelse(k <= 0, NA, C0)
  )
fallback_k <- min(start_vals_lin_clean$k_clean, na.rm = TRUE)

C0_start_by_soil <- start_vals_lin_clean %>%
  filter(k > 0) %>%
  group_by(Locality) %>%
  summarize(C0 = median(C0, na.rm = TRUE)) %>%
  pull(C0)

C0_start   <- median(C0_start_by_soil, na.rm = TRUE)
logk_start <- log(median(start_vals_lin_clean$k_clean, na.rm = TRUE))

# ---- nested non-linear mixed effect models ----
#1. Single temperature sensitivity, single initial available organic carbon concentration
M0 <- nlme(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = Ctr,
  fixed = list(C0 ~ 1, logk ~ Temp),
  random = C0 ~ 1 | Locality,
  method = "ML",
  start = c(C0 = C0_start, logk = c(logk_start, 0)),
  control = nlmeControl(maxIter = 300, msMaxIter = 300, msMaxEval = 500),
  na.action = na.exclude
)
summary(M0)
#2. Initial available organic carbon concentration differs among soils - almost certainty
M1 <- update(M0, fixed = list(C0 ~ Locality, logk ~ Temp),
             start = c(C0 = rep(1078.8, 6), logk = c(-5.6469, 0.1808)))
anova(M0, M1)
summary(M1)

# ---- random variability component dropped to almost zero so I will refit without random component and compare ----
M1_gnls <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = Ctr,
  params = list(C0 ~ Locality, logk ~ Temp),
  start = c(C0 = c(1824.16, -76.99, -1019.85, -791.49, -1376.17, -946.77),
            logk = c(-5.7250, 0.1758)),
  na.action = na.exclude
)

summary(M1_gnls)
AIC(M1, M1_gnls)

# ---- fitting heteroscedascity on M1 ----
M1_diag <- tibble(fitted = fitted(M1_gnls), resid = residuals(M1_gnls),
                  Temp = Ctr$Temp, Locality = Ctr$Locality)

M1_varIdent_temp <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = Ctr,
  params = list(C0 ~ Locality, logk ~ Temp),
  start = coef(M1_gnls),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300),
  na.action = na.exclude
)
anova(M1_gnls, M1_varIdent_temp)
AIC(M1_gnls, M1_varIdent_temp)

anova(M1_varIdent_temp)
R2M(M1_varIdent_temp)

# ======================= Non-linear model fitting with initial SOC defined  - nested models 
toc_anchor_soil <- Ctr %>% group_by(Locality) %>%
  summarize(TOC_anchor = max(soil_C*1000, na.rm=TRUE))   # pooled across all temp/FT, same units as SoilResp

Ctr_cc <- Ctr[complete.cases(Ctr[, c("SoilResp", "Temp", "Time", "Locality")]), ]
nrow(Ctr_cc)

Ctr_cc <- Ctr_cc %>% left_join(toc_anchor_soil, by = "Locality")

M1_aB <- gnls(
  SoilResp ~ exp(logk) * plogis(logit_a) * TOC_anchor * exp(-exp(logk) * Time),
  data = Ctr_cc,
  params = list(logk ~ Temp, logit_a ~ 1),
  start = c(logk = coef(M1_gnls)[c("logk.(Intercept)","logk.Temp")], logit_a = qlogis(0.01)),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)


# ---- This is now the null model with a fraction of total organic C directly specified in the equation ----
N0 <- gnls(
  SoilResp ~ exp(logk) * plogis(logit_a) * TOC_anchor * exp(-exp(logk) * Time),
  data = Ctr_cc,
  params = list(logk ~ Temp, logit_a ~ 1),
  start = c(logk = coef(M1_aB)[c("logk.(Intercept)","logk.Temp")], logit_a = coef(M1_aB)["logit_a"]),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)

summary(N0)

# ---- Does the locality affect a? According to all we know, it should ----
N1 <- update(N0, params = list(logk ~ Temp, logit_a ~ Locality),
             start = c(coef(N0)[1:2], rep(coef(N0)["logit_a"], 6)))

summary(N1)
anova(N0, N1)
AIC(N0, N1)

# ---- Does the temperature affect a? ----
N2 <- update(N1, params = list(logk ~ Temp, logit_a ~ Locality + Temp),
             start = c(coef(N1), 0))

summary(N2)
anova(N1, N2)
AIC(N1, N2)
R2M(N2)
anova(N2, type = 'sequential', adjustSigma = T)

# ---- The next model is already complex, it tests the interaction between the locality and temperature ----
N3 <- update(N2, params = list(logk ~ Temp, logit_a ~ Locality * Temp),
             start = c(coef(N2)[1:8], rep(0, 6)))

summary(N3)
anova(N2, N3)
AIC(N2, N3)

# ---- Is the kinetic constant locality-specific? ----
N2k1 <- update(N2, params = list(logk ~ Locality + Temp, logit_a ~ Locality + Temp),
               start = c(rep(coef(N2)["logk.(Intercept)"], 6), coef(N2)["logk.Temp"],
                         coef(N2)[3:9]))

summary(N2k1)
anova(N2, N2k1)
AIC(N2, N2k1)

#===============================================================================
# ------------------- Results of the best model N2----------------------------->
#===============================================================================

# ---- Bootstraping the best model N2 ----
set.seed(123)
N2_boot <- boot_gnls(N2, nboot = 999)
boot_coefs_N2 <- as.data.frame(N2_boot$t)
names(boot_coefs_N2) <- names(N2_boot$t0)
boot_coefs_N2 <- boot_coefs_N2[complete.cases(boot_coefs_N2), ]
nrow(boot_coefs_N2)

# ---- calculating Q10 value of k ----
logk_Temp_draws <- boot_coefs_N2[["logk.Temp"]]
Q10_soil_final <- exp(10 * logk_Temp_draws)

tibble(
  metric = "Q10 (soil-derived, availability-corrected)",
  mean = mean(Q10_soil_final), median = median(Q10_soil_final),
  lwr = quantile(Q10_soil_final, 0.025), upr = quantile(Q10_soil_final, 0.975)
)

# ---- available fraction of TOC ----
localities <- levels(Ctr_cc$Locality)
temp_grid <- seq(-5, 5, by=1)

a_traj <- lapply(localities, function(loc) {
  base <- boot_coefs_N2[["logit_a.(Intercept)"]]
  loc_off <- if (loc == localities[1]) 0 else boot_coefs_N2[[paste0("logit_a.Locality", loc)]]
  temp_coef <- boot_coefs_N2[["logit_a.Temp"]]
  
  sapply(temp_grid, function(T) {
    a_draws <- plogis(base + loc_off + temp_coef * T)
    c(mean = median(a_draws), lwr = quantile(a_draws,0.025), upr = quantile(a_draws,0.975))
  }) %>% t() %>% as_tibble() %>% mutate(Temp = temp_grid, Locality = loc)
}) %>% bind_rows()

print(a_traj)

# ---- For different localities ----
a_loc <- lapply(localities, function(loc) {
  base <- boot_coefs_N2[["logit_a.(Intercept)"]]
  loc_off <- if (loc == localities[1]) 0 else boot_coefs_N2[[paste0("logit_a.Locality", loc)]]
  temp_coef <- boot_coefs_N2[["logit_a.Temp"]]
  a_draws <- plogis(base + loc_off + + temp_coef * 0)*100
  c(mean = median(a_draws), lwr = quantile(a_draws,0.025), upr = quantile(a_draws,0.975)) %>%
    t() %>% as_tibble() %>% mutate(Locality = loc)
}) %>% bind_rows()

print(a_loc) 

# ---- fold change ratio between -5 and +5 ----
a_neg5 <- a_N2fun(boot_coefs_N2, "RA", -5)
a_5    <- a_N2fun(boot_coefs_N2, "RA", 5)

a_ratio <- a_5 / a_neg5

tibble(
  metric = "Fold-change in available fraction, -5 to 5C",
  median = median(a_ratio),
  lwr = quantile(a_ratio, 0.025),
  upr = quantile(a_ratio, 0.975)
)

# ---- Carbon loss estimate ----
RootFreeLoss <- data.frame(Locality = rep(unique(TrLongerFull$Locality), each = 3),
                           Temp = rep(c(-5, 0, 5), times = 6),
                           Center = numeric(length = 6*3),
                           lw = numeric(length = 6*3),
                           up = numeric(length = 6*3))
RootFreeLoss$Locality <- factor(RootFreeLoss$Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA"))
for(i in unique(RootFreeLoss$Locality)){
  for(j in c(-5, 0, 5)){
    closs <- a_N2fun(boot_coefs_N2, i, j)*(1-exp(-logk_N2fun(boot_coefs_N2, j)*120))*100
    RootFreeLoss[(RootFreeLoss$Locality == i & RootFreeLoss$Temp == j), 3:5] <- c(mean(closs),
                                                                                  quantile(closs, probs = c(0.025, 0.975)))
  }
}
for(i in  c("SC", "VW", "NW", "FU", "JA", "RA")){
  print(RootFreeLoss[RootFreeLoss$Locality == i, ])
}

#===============================================================================
# 2.Freeze-thaw cycles (-5/5, -18/5) + soil respiration only
# - the best model N2 is used to make predictions on respiration in 
#   alternating thermal regimes after 2 months of incubation
# - metanalytical approach is used to quantify difference between predictions and 
#   measurements
# - the trajectory of freeze-thaw effect is extrapolated to later incubations
#===============================================================================
# ---- We start with -5/5 treatment because it is within range of stable temperatures that N2 is calibrated against (no extrapolation to -18) ----
# ---- Model parameters
logk_N2fun <- function(draws, temp) {
  exp(draws[["logk.(Intercept)"]] + draws[["logk.Temp"]] * temp)
}

a_N2fun <- function(draws, locality, temp) {
  base <- draws[["logit_a.(Intercept)"]]
  off  <- if (locality == localities[1]) 0 else draws[[paste0("logit_a.Locality", locality)]]
  plogis(base + off + draws[["logit_a.Temp"]] * temp)
}

# ---- We will set the reference point at month 2 for -5/5 treatment. The accuracy of the reference point is cross-validated against the same predictions performed for stable -5 and +5 ---- 
k_neg5 <- logk_N2fun(boot_coefs_N2, -5)
k_5    <- logk_N2fun(boot_coefs_N2, 5)

predict_month2 <- function(locality) {
  toc_val <- toc_anchor_soil$TOC_anchor[toc_anchor_soil$Locality == locality]
  a_neg5 <- a_N2fun(boot_coefs_N2, locality, -5)
  a_5    <- a_N2fun(boot_coefs_N2, locality, 5)
  
  pred_stable_neg5 <- k_neg5 * a_neg5 * toc_val * exp(-k_neg5 * 60.88)
  pred_stable_5    <- k_5    * a_5    * toc_val * exp(-k_5 * 60.88)
  pred_5_5         <- k_5    * a_5    * toc_val * exp(-(k_neg5 * 60.88 + k_5 * 14))
  
  tibble(
    Series = c("stable_-5","stable_5","-5/5"),
    pred_mean = c(mean(pred_stable_neg5), mean(pred_stable_5), mean(pred_5_5)),
    pred_lwr  = c(quantile(pred_stable_neg5,0.025), quantile(pred_stable_5,0.025), quantile(pred_5_5,0.025)),
    pred_upr  = c(quantile(pred_stable_neg5,0.975), quantile(pred_stable_5,0.975), quantile(pred_5_5,0.975)),
    Locality = locality
  )
}

pred_month2_all <- lapply(localities, predict_month2) %>% bind_rows()

# ---- Comparing with observations ----
FT_df <- TrLongerFull %>%
  filter(Origin == 'Control soil' & TempFactor %in% c("5", "-18/5", "-5/5")) %>%
  mutate(TempFactor = factor(TempFactor, levels = c("5", "-18/5", "-5/5")),
         Locality = factor(Locality))
FT_df$Time <- FT_df$MonthsInc*30

obs_FT <- FT_df %>%
  filter(TempFactor %in% c("-5/5", "-18/5")) %>%
  mutate(Time_since_thaw = round((Time - 60.88) / 30.44)) %>%   # months since thaw, from days
  rename(Treatment = TempFactor)

obs_month2 <- bind_rows(
  Ctr_cc %>% filter(Temp==-5, round((Time-60.88)/30.44)==0) %>% mutate(Series="stable_-5"),
  Ctr_cc %>% filter(Temp==5,  round((Time-60.88)/30.44)==0) %>% mutate(Series="stable_5"),
  obs_FT %>% filter(Treatment=="-5/5", Time_since_thaw==0) %>% mutate(Series="-5/5")
) %>%
  group_by(Series, Locality) %>%
  summarize(obs_mean = mean(SoilResp, na.rm=TRUE), .groups="drop")

compare_month2 <- left_join(obs_month2, pred_month2_all, by=c("Series","Locality"))
print(compare_month2, n=30)

# ---- Formal test ----
compare_month2 <- compare_month2 %>% mutate(within_band = obs_mean >= pred_lwr & obs_mean <= pred_upr)
compare_month2 %>% group_by(Series) %>% summarize(n_within = sum(within_band), n_total = n())

# ---- accounting for variance in measurements and predictions ----
# Observations
obs_summary <- obs_FT %>% filter(Treatment == "-5/5", Time_since_thaw == 0) %>%
  group_by(Locality) %>%
  summarize(obs_mean = mean(SoilResp, na.rm=TRUE),
            obs_sd = sd(SoilResp, na.rm=TRUE),
            n = n(), .groups = "drop")

# Calculating an excess with the proper variance 
meta_data <- pred_5_5_val %>%   # control-soil predicted values, from the earlier N2-based check
  filter(Root == "absent") %>%
  left_join(obs_summary, by = "Locality") %>%
  mutate(
    obs_se = obs_sd / sqrt(n),
    pred_se = (pred_upr - pred_lwr) / (2*1.96),
    yi = log(obs_mean) - log(pred_mean),
    vi = (obs_se/obs_mean)^2 + (pred_se/pred_mean)^2
  )

print(meta_data)

# Formal test
res_overall <- rma(yi, vi, data = meta_data, method = "REML")
summary(res_overall)

# back-transformed pooled excess ratio
exp(c(estimate = coef(res_overall), ci.lb = res_overall$ci.lb, ci.ub = res_overall$ci.ub))

# try elevation modifier
elev_lookup <- tibble(Locality = c("SC","VW","NW","FU","RA","JA"),
                      Elevation = c(2160,2160,2130,1780,1470,1330),
                      Vegetation = c('Scree', "Meadow", "Meadow", "Meadow", 'Forest', 'Forest'))

meta_data_elev <- meta_data %>% left_join(elev_lookup, by="Locality")
res_elev <- rma(yi, vi, mods = ~ Elevation, data = meta_data_elev, method="REML")
summary(res_elev)
res_veg <- rma(yi, vi, mods = ~ Vegetation, data = meta_data_elev, method="REML")
summary(res_veg)

# ---- the same for -18/5 ----
# ---- accounting for variance in measurements and predictions ----
# Observations
obs_summary_18 <- obs_FT %>% filter(Treatment == "-18/5", Time_since_thaw == 0) %>%
  group_by(Locality) %>%
  summarize(obs_mean = mean(SoilResp, na.rm=TRUE),
            obs_sd = sd(SoilResp, na.rm=TRUE),
            n = n(), .groups = "drop")

# Calculating an excess with the proper variance 
meta_data_18 <- pred_month2_FT %>%   # control-soil predicted values, from the earlier N2-based check
  filter(TempFactor == "-18/5") %>%
  left_join(obs_summary_18, by = "Locality") %>%
  mutate(
    obs_se = obs_sd / sqrt(n),
    pred_se = (pred_upr - pred_lwr) / (2*1.96),
    yi = log(obs_mean) - log(pred_mean),
    vi = (obs_se/obs_mean)^2 + (pred_se/pred_mean)^2
  )

print(meta_data_18)

# Formal test
res_overall_18 <- rma(yi, vi, data = meta_data_18, method = "REML")
summary(res_overall_18)

# back-transformed pooled excess ratio
exp(c(estimate = coef(res_overall_18), ci.lb = res_overall_18$ci.lb, ci.ub = res_overall_18$ci.ub))

# try elevation modifier
meta_data_elev_18 <- meta_data_18 %>% left_join(elev_lookup, by="Locality")
res_elev_18 <- rma(yi, vi, mods = ~ Elevation, data = meta_data_elev_18, method="REML")
summary(res_elev_18)
res_veg_18 <- rma(yi, vi, mods = ~ Vegetation, data = meta_data_elev_18, method="REML")
summary(res_veg_18)

# ---- Now, we will compare time trajectories against the well established reference point ----
rel_decline <- bind_rows(
  obs_FT %>% filter(Treatment == "-5/5") %>% mutate(Series = "-5/5"),
  Ctr_cc %>% filter(Temp == -5) %>% transmute(Locality, Time_since_thaw = round((Time-60.88)/30.44), SoilResp, Series = "stable_-5"),
  Ctr_cc %>% filter(Temp == 5)  %>% transmute(Locality, Time_since_thaw = round((Time-60.88)/30.44), SoilResp, Series = "stable_5")
) %>%
  group_by(Series, Locality) %>%
  mutate(month2_val = mean(SoilResp[Time_since_thaw == 0], na.rm = TRUE)) %>%
  filter(!is.na(month2_val), month2_val > 0) %>%
  mutate(rel_change = SoilResp / month2_val) %>%
  ungroup()

rel_decline %>% count(Series, Time_since_thaw)  

rel_decline %>% group_by(Series, Time_since_thaw) %>%
  summarize(mean_rel = mean(rel_change, na.rm=TRUE), se = sd(rel_change,na.rm=TRUE)/sqrt(n()), .groups="drop")

#==============time to excess respiration attenuation
set.seed(123)
n_boot <- 2000

# 1. Bootstrap draws of the month-0 excess ratio, from the rma() model's own sampling distribution
ratio0_draws <- exp(rnorm(n_boot, mean = coef(res_overall), sd = res_overall$se))

# 2. Bootstrap the relative-decline trajectory by resampling replicates within each Series/Time cell
rel_boot <- replicate(n_boot, {
  boot_sample <- rel_decline %>% filter(Series == "-5/5") %>%
    group_by(Locality, Time_since_thaw) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    summarize(mean_rel = mean(rel_change, na.rm=TRUE), .groups="drop")
  boot_sample$mean_rel
}) %>% t()   # n_boot x 3 matrix, columns = Time_since_thaw 0,1,2

colnames(rel_boot) <- unlist(lapply(unique(TrLongerFull$Locality), paste0, c("-m0", "-m2", "-m3")))


# 3. Combine: ratio(m) = ratio(0) * rel_obs(m)/rel_obs(0), but rel_obs(0) = 1 by construction
ratio_traj_boot <- ratio0_draws * rel_boot   # element-wise: each row is one bootstrap replicate's full trajectory

ratio_traj_summary <- tibble(
  Legend = colnames(rel_boot),
  median = apply(ratio_traj_boot, 2, median, na.rm = T),
  lwr = apply(ratio_traj_boot, 2, quantile, 0.025, na.rm = T),
  upr = apply(ratio_traj_boot, 2, quantile, 0.975, na.rm = T)
)
print(ratio_traj_summary)

# boot_coefs_N2 and ratio0_draws/rel_boot have different lengths -> resample to align
# (treats the two models' bootstraps as independent)
set.seed(123)
idx <- sample(nrow(boot_coefs_N2), n_boot, replace = TRUE)
kn <- logk_N2fun(boot_coefs_N2, -5)[idx]
kp <- logk_N2fun(boot_coefs_N2,  5)[idx]

scen <- list(all_cold = kn * 30,
             all_warm = kp * 30)

traj <- lapply(scen, function(e) ratio0_draws * rel_boot / cbind(1, exp(-e), exp(-2*e),
                                                                 1, exp(-e), exp(-2*e),
                                                                 1, exp(-e), exp(-2*e),
                                                                 1, exp(-e), exp(-2*e),
                                                                 1, exp(-e), exp(-2*e),
                                                                 1, exp(-e), exp(-2*e)))

trajDFrame <- data.frame(Scenario = rep(c("Cold", "Warm"), each = 2000*18),
                         Locality = rep(rep(unique(TrLongerFull$Locality), each = 3*2000), 2),
                         Month = rep(rep(c(0, 1, 2), each = 2000), times = 6*2),
                         Rm = as.numeric(c(traj$all_cold, traj$all_warm))) 

trajDFrame |> group_by(Scenario, Month) |> summarise(Center = mean(Rm, na.rm = T),
                                                     lb = quantile(Rm, probs = c(0.025), na.rm = T),
                                                     ub = quantile(Rm, probs = c(0.975), na.rm = T))

# ---- adding -18/5 treatment ----
rel_decline_full <- bind_rows(
  rel_decline,   # already has -5/5, stable_-5, stable_5
  obs_FT %>% filter(Treatment == "-18/5") %>% mutate(Series = "-18/5") %>%
    group_by(Locality) %>%
    mutate(month2_val = mean(SoilResp[Time_since_thaw == 0], na.rm = TRUE)) %>%
    filter(!is.na(month2_val), month2_val > 0) %>%
    mutate(rel_change = SoilResp / month2_val) %>%
    ungroup()
)

rel_decline_full %>% count(Series, Time_since_thaw)

rel_decline_full %>% filter(Series %in% c("-5/5","-18/5")) %>%
  group_by(Locality, Time_since_thaw) %>%
  summarize(mean_rel = mean(rel_change, na.rm=TRUE), .groups="drop") %>%
  tidyr::pivot_wider(names_from = Time_since_thaw, values_from = mean_rel)
# 
# rel_model_full <- lm(rel_change ~ Series * Time_since_thaw, data = rel_decline_full %>% filter(Time_since_thaw > 0))
# anova(rel_model_full)
# summary(rel_model_full)

#==============time to excess respiration attenuation
set.seed(123)
n_boot <- 2000

# 1. Bootstrap draws of the month-0 excess ratio, from the rma() model's own sampling distribution
ratio0_draws_18 <- exp(rnorm(n_boot, mean = coef(res_overall_18), sd = res_overall_18$se))

# 2. Bootstrap the relative-decline trajectory by resampling replicates within each Series/Time cell
rel_boot_18 <- replicate(n_boot, {
  boot_sample <- rel_decline_full %>% filter(Series == "-18/5") %>%
    group_by(Locality, Time_since_thaw) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    summarize(mean_rel = mean(rel_change, na.rm=TRUE), .groups="drop")
  boot_sample$mean_rel
}) %>% t()   # n_boot x 3 matrix, columns = Time_since_thaw 0,1,2

colnames(rel_boot_18) <- unlist(lapply(unique(TrLongerFull$Locality), paste0, c("-m0", "-m2", "-m3")))

# 3. Combine: ratio(m) = ratio(0) * rel_obs(m)/rel_obs(0), but rel_obs(0) = 1 by construction
ratio_traj_boot_18 <- ratio0_draws_18 * rel_boot_18   # element-wise: each row is one bootstrap replicate's full trajectory

ratio_traj_summary_18 <- tibble(
  Legend = colnames(rel_boot_18),
  median = apply(ratio_traj_boot_18, 2, median, na.rm = T),
  lwr = apply(ratio_traj_boot_18, 2, quantile, 0.025, na.rm = T),
  upr = apply(ratio_traj_boot_18, 2, quantile, 0.975, na.rm = T)
)
print(ratio_traj_summary_18)

# boot_coefs_N2 and ratio0_draws/rel_boot have different lengths -> resample to align
# (treats the two models' bootstraps as independent)
set.seed(123)
idx <- sample(nrow(boot_coefs_N2), n_boot, replace = TRUE)
kn <- logk_N2fun(boot_coefs_N2, -5)[idx]
kp <- logk_N2fun(boot_coefs_N2,  5)[idx]

scen18 <- list(all_cold = kn * 0,
               all_warm = kp * 30)

traj18 <- lapply(scen18, function(e) ratio0_draws_18 * rel_boot_18 / cbind(1, exp(-e), exp(-2*e),
                                                                           1, exp(-e), exp(-2*e),
                                                                           1, exp(-e), exp(-2*e),
                                                                           1, exp(-e), exp(-2*e),
                                                                           1, exp(-e), exp(-2*e),
                                                                           1, exp(-e), exp(-2*e)))

trajDFrame18 <- data.frame(Scenario = rep(c("Cold", "Warm"), each = 2000*18),
                           Locality = rep(rep(unique(TrLongerFull$Locality), each = 3*2000), 2),
                           Month = rep(rep(c(0, 1, 2), each = 2000), times = 6*2),
                           Rm = as.numeric(c(traj18$all_cold, traj18$all_warm))) 

trajDFrame18 |> group_by(Scenario, Month) |> summarise(Center = mean(Rm, na.rm = T),
                                                       lb = quantile(Rm, probs = c(0.025), na.rm = T),
                                                       ub = quantile(Rm, probs = c(0.975), na.rm = T))


# ---- Direct comparison of measurements at month 4 ----
#1. -5/5
anova(aov(SoilResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Roots == 'No')))
anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Roots == 'No')))
summary(aov(SoilResp ~ TempFactor, data = TrLongerFull,
            subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Roots == 'No')))
TukeyHSD(aov(SoilResp ~ TempFactor, data = TrLongerFull,
             subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Roots == 'No')))

#1. -18/5
anova(aov(SoilResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Roots == 'No')))
anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Roots == 'No')))
summary(aov(SoilResp ~ TempFactor, data = TrLongerFull,
            subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Roots == 'No')))
TukeyHSD(aov(SoilResp ~ TempFactor, data = TrLongerFull,
             subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Roots == 'No')))
TukeyHSD(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
             subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Roots == 'No')))

Soil18end <- as.data.frame(TukeyHSD(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
                                        subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5' & Roots == 'No')))$`TempFactor:Locality`)

for(i in unique(TrLonger$Locality)){
  #Soil specific comparison
  comp5negative <- paste0('-5:', i, '--18/5:', i)
  comp5positive <- paste0('5:', i, '--18/5:', i)
  
  print(rbind(
    Soil18end[grepl(comp5negative, Soil18end$comp), c("diff", 'p adj', 'comp')],
    Soil18end[grepl(comp5positive, Soil18end$comp), c("diff", 'p adj', 'comp')]
  )[2:3, ])
  
}

#===============================================================================
# Figures for root-free soils
#===============================================================================
#=======================Stable temperatures - measured and predicted 
#Measured
RootFreeStable <- Ctr_cc |> group_by(Locality, Temp, Time) |> summarise(y = mean(SoilResp, na.rm = T),
                                                                        ySD = sd(SoilResp, na.rm = T),
                                                                        TOC_anchor = median(TOC_anchor, na.rm = T))
#Predicted
RootFreeStable$yhat <- numeric(length = nrow(RootFreeStable))
RootFreeStable$lwr <- numeric(length = nrow(RootFreeStable))
RootFreeStable$upr <- numeric(length = nrow(RootFreeStable))

for(i in 1:nrow(RootFreeStable)){
  #Model parameters for row i
  loc_i <- RootFreeStable$Locality[i]
  k_i = exp(boot_coefs_N2[['logk.(Intercept)']] + boot_coefs_N2[['logk.Temp']]*RootFreeStable$Temp[i])
  if(loc_i == 'FU'){
    alfa_i =  plogis(boot_coefs_N2[['logit_a.(Intercept)']] + 
                       boot_coefs_N2[['logit_a.Temp']]*RootFreeStable$Temp[i])
  } else {
    alfa_i = plogis(boot_coefs_N2[['logit_a.(Intercept)']] + 
                      boot_coefs_N2[[paste0('logit_a.Locality', loc_i)]] + 
                      boot_coefs_N2[['logit_a.Temp']]*RootFreeStable$Temp[i])
  }
  
  #Prediction for row i
  RootFreeStable[i, c('yhat', 'lwr', 'upr')] <- matrix(c(mean(RootFreeStable$TOC_anchor[i]*
                                                                k_i*alfa_i*exp(-k_i*RootFreeStable$Time[i])), 
                                                         quantile(x = RootFreeStable$TOC_anchor[i]*
                                                                    k_i*alfa_i*exp(-k_i*RootFreeStable$Time[i]),
                                                                  probs = c(0.025, 0.975))), nrow = 1)
  
  
}

RootFreeStable_Figure <- RootFreeStable %>%
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
       y = expression(paste("Respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
  new_scale_fill() +
  new_scale_color() +
  geom_point(cex = 4, pch = 21, aes(y = yhat, fill = factor(Temp)), position = position_dodge(width = 1),
             show.legend = F) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr, colour = factor(Temp)), 
                position = position_dodge(width = 1), width = 0.1, show.legend = F) +
  scale_fill_manual(values = c('red', 'red', 'red')) +
  scale_color_manual(values = c('red', 'red', 'red')) +
  theme(legend.position = c(0.08, 0.88))

ggsave("StatOut/RootFreeStable_Figure.pdf", RootFreeStable_Figure,
  width = 300, height = 180, units = "mm", device = cairo_pdf)

#=======================Freeze-thaw cycles - temporal trajectories 
trajDFrame$Series <- c('-5/5')
trajDFrame18$Series <- c('-18/5')

RootFreeFT_temporalFig <- rbind(trajDFrame, trajDFrame18) |> group_by(Series, Scenario, Locality, Month) |> 
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

ggsave("StatOut/RootFreeFT_temporalFig.pdf", RootFreeFT_temporalFig,
       width = 360, height = 180, units = "mm", device = cairo_pdf)

#=======================Freeze-thaw cycles - composite figure
predict_month2_FT <- function(locality) {
  toc_val <- toc_anchor_soil$TOC_anchor[toc_anchor_soil$Locality == locality]
  a_neg5 <- a_N2fun(boot_coefs_N2, locality, -5)
  a_5    <- a_N2fun(boot_coefs_N2, locality, 5)
  
  pred_18_5        <- k_5    * a_5    * toc_val * exp(-(k_5 * 14))
  pred_5_5         <- k_5    * a_5    * toc_val * exp(-(k_neg5 * 60.88 + k_5 * 14))
  
  tibble(
    TempFactor = c("-18/5","-5/5"),
    pred_mean = c(mean(pred_18_5), median(pred_5_5)),
    pred_lwr  = c(quantile(pred_18_5,0.025), quantile(pred_5_5,0.025)),
    pred_upr  = c(quantile(pred_18_5,0.975), quantile(pred_5_5,0.975)),
    Locality = locality
  )
}

pred_month2_FT <- lapply(localities, predict_month2_FT) %>% bind_rows()

RootFreeRef <- Ctr_cc |> group_by(Locality, Temp, MonthsInc) |> 
  summarise(yref = mean(SoilResp, na.rm = T), yrefSD = sd(SoilResp, na.rm = T)) |> 
  filter(Temp != 0) |> filter(MonthsInc != 3) |>  
  mutate(TempF = factor(Temp, levels = c(-5, 5), labels = c("-5˚C", '+5˚C')))

pred_month2_FT2 <- rbind(pred_month2_FT, pred_month2_FT)
pred_month2_FT2$MonthsInc <- rep(c(2, 4), each = nrow(pred_month2_FT)) 
pred_month2_FT2[(nrow(pred_month2_FT)+1):nrow(pred_month2_FT2), c(2:4)] <- NA

RootFreeFT_compositeFig <- TrLongerFull |> filter(Roots == 'No' & StabTemp == 'No' & MonthsInc != 3) |> 
  group_by(Locality, TempFactor, MonthsInc) |> summarise(y = mean(SoilResp, na.rm = T),
                                                         ySD = sd(SoilResp, na.rm = T)) |> 
  inner_join(pred_month2_FT2, by = c('Locality', "TempFactor", "MonthsInc")) |> 
  full_join(RootFreeRef) |> 
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
       y = expression(paste("Respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
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

ggsave("StatOut/RootFreeFT_compositeFig.pdf", RootFreeFT_compositeFig,
       width = 320, height = 180, units = "mm", device = cairo_pdf)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF SCRIPT