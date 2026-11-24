#===============================================================================
#                               Tatra Mountains 
#                             Incubation experiment
#
# This script analyzes root amended soils:
#       - root-derived respiration (SOC-derived respiration is in different file)
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
# 1. Stable temperatures (-5, 0, 5) + root-derived respiration only [non-linear mixed effects]
# - nested models are designed and compared 
# - the best one is selected
#===============================================================================
# ======================= Data preparation 
TrRoots <- TrLongerFull |> filter(Origin == "Roots" & StabTemp == 'Yes')
TrRoots$Locality <- as.factor(TrRoots$Locality)
# ---- Converting C0 from mg to ug as MonthInc from months to Time in days
TrRoots$C0 <- TrRoots$C0*1000 #ug per gram of dry soil
TrRoots$Time <- TrRoots$MonthsInc*30 #days
#checking the ranges
cbind(summary(TrRoots$RootResp), summary(TrRoots$C0))
TrRoots_cc <- TrRoots[complete.cases(TrRoots[, c("SoilResp", "Temp", "Time", "Locality")]), ]
nrow(TrRoots_cc)

# ======================= Non-linear model fitting - nested models 
# ---- the null model with available fraction of root material directly acknowledged ---- 
Pa_k0 <- gnls(
  RootResp ~ exp(logk) * plogis(logit_a) * C0 * exp(-exp(logk) * Time),
  data = TrRoots_cc,
  params = list(logk ~ 1, logit_a ~ 1),
  start = c(logk = log(0.05), logit_a = qlogis(0.15)),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)
# ---- Exponential temperature function ----
Pa_k1 <- update(Pa_k0, params = list(logk ~ Temp, logit_a ~ 1),
                start = c(coef(Pa_k0)["logk"], 0, coef(Pa_k0)["logit_a"]))
# ---- Gaussian temperature function ----
Pa_k2 <- update(Pa_k1, params = list(logk ~ Temp + I(Temp^2), logit_a ~ 1),
                start = c(coef(Pa_k1)[1:2], 0, coef(Pa_k1)[3]))

anova(Pa_k0, Pa_k1, Pa_k2)
AIC(Pa_k0, Pa_k1, Pa_k2)
summary(Pa_k2)

# ---- the statistical preference of Gaussian temperature response function is tight. Therefore, logit_a is first allowed to explain some curvature before decding for Gaussian ----
# ---- with exponential function ----
Pb_k1_aT <- gnls(
  RootResp ~ exp(logk) * plogis(logit_a) * C0 * exp(-exp(logk) * Time),
  data = TrRoots_cc,
  params = list(logk ~ Temp, logit_a ~ Temp),
  start = c(logk = coef(Pa_k1)[1:2], logit_a = c(coef(Pa_k1)[3], 0)),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)

summary(Pb_k1_aT)

# ---- with the Gaussian function ----
Pb_k2_aT <- update(Pb_k1_aT, params = list(logk ~ Temp + I(Temp^2), logit_a ~ Temp),
                   start = c(coef(Pb_k1_aT)[1:2], 0, coef(Pb_k1_aT)[3:4]))

anova(Pb_k1_aT, Pb_k2_aT)
AIC(Pb_k1_aT, Pb_k2_aT)

AIC(Pa_k0, Pa_k1, Pa_k2, Pb_k1_aT, Pb_k2_aT)

# ---- does available fraction of root material changes with the locality? ----
Pc_a_loc <- update(Pb_k1_aT, params = list(logk ~ Temp, logit_a ~ Locality + Temp),
                   start = c(coef(Pb_k1_aT)[1:2], coef(Pb_k1_aT)["logit_a.(Intercept)"],
                             rep(0,5), coef(Pb_k1_aT)["logit_a.Temp"]))

anova(Pb_k1_aT, Pc_a_loc)
AIC(Pb_k1_aT, Pc_a_loc)

# ---- does k of root material changes with the locality? ----
Pc_k_loc <- update(Pb_k1_aT, params = list(logk ~ Temp + Locality, logit_a ~ Temp),
                   start = c(coef(Pb_k1_aT)[1], rep(0,5), coef(Pb_k1_aT)[2], coef(Pb_k1_aT)["logit_a.(Intercept)"],
                             coef(Pb_k1_aT)["logit_a.Temp"]))

anova(Pb_k1_aT, Pc_k_loc) 
AIC(Pb_k1_aT, Pc_k_loc)

AIC(Pb_k1_aT, Pc_a_loc, Pc_k_loc)

#===============================================================================
# ----------------- Results of the best model Pc_k_loc------------------------->
#===============================================================================
R2M(Pc_k_loc)
anova(Pc_k_loc, type = 'marginal')

Pc_k_loc_boot <- boot_gnls(Pc_k_loc, nboot=999)
boot_coefs_Pc_k <- as.data.frame(Pc_k_loc_boot$t) 
names(boot_coefs_Pc_k) <- names(Pc_k_loc_boot$t0)
boot_coefs_Pc_k <- boot_coefs_Pc_k[complete.cases(boot_coefs_Pc_k),]

# ---- Q10 ----
logk_Temp_draws <- boot_coefs_Pc_k[["logk.Temp"]]
Q10_root_final <- exp(10 * logk_Temp_draws)

tibble(
  metric = "Q10 (root-derived, availability-corrected)",
  median = median(Q10_root_final),
  lwr = quantile(Q10_root_final, 0.025),
  upr = quantile(Q10_root_final, 0.975)
)

# ---- k(T), locality-specific ----
k_traj <- lapply(localities, function(loc) {
  sapply(temp_grid, function(T) {
    k_draws <- k_Pck_fun(boot_coefs_Pc_k, loc, T)
    c(median = median(k_draws), lwr = quantile(k_draws, 0.025), upr = quantile(k_draws, 0.975))
  }) %>% t() %>% as_tibble() %>% mutate(Temp = temp_grid, Locality = loc)
}) %>% bind_rows()

print(k_traj) |> filter(Temp == 0)

# ---- available fraction of TOC ----
aRoot_traj <- sapply(temp_grid, function(T) {
  a_draws <- a_Pck_fun(boot_coefs_Pc_k, T)
  c(median = median(a_draws), lwr = quantile(a_draws, 0.025), upr = quantile(a_draws, 0.975))
}) %>% t() %>% as_tibble() %>% mutate(Temp = temp_grid)

print(aRoot_traj)

# ---- fold change ratio between -5 and +5 ----
aRoots_neg5 <- a_Pck_fun(boot_coefs_Pc_k, -5)
aRoots_5    <- a_Pck_fun(boot_coefs_Pc_k, 5)

aRoots_ratio <- aRoots_5 / aRoots_neg5

tibble(
  metric = "Fold-change in available fraction, -5 to 5C",
  median = median(aRoots_ratio),
  lwr = quantile(aRoots_ratio, 0.025),
  upr = quantile(aRoots_ratio, 0.975)
)

# ---- Carbon loss estimate ----
RootCLoss <- data.frame(Locality = rep(unique(TrLongerFull$Locality), each = 3),
                        Temp = rep(c(-5, 0, 5), times = 6),
                        Center = numeric(length = 6*3),
                        lw = numeric(length = 6*3),
                        up = numeric(length = 6*3))
RootCLoss$Locality <- factor(RootCLoss$Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA"))

for(i in unique(RootCLoss$Locality)){
  for(j in c(-5, 0, 5)){
    closs <- a_Pck_fun(boot_coefs_Pc_k, j)*(1-exp(-k_Pck_fun(boot_coefs_Pc_k, i, j)*120))*100
    RootCLoss[(RootCLoss$Locality == i & RootCLoss$Temp == j), 3:5] <- c(median(closs),
                                                                         quantile(closs, probs = c(0.025, 0.975)))
  }
}
for(i in  c("SC", "VW", "NW", "FU", "JA", "RA")){
  print(RootCLoss[RootCLoss$Locality == i, ])
}

#===============================================================================
# 2.Freeze-thaw cycles (-5/5, -18/5) + Root-derived respiration in root-amended soils only
# - the best model Pc_k_loc is used to make predictions on respiration in 
#   alternating thermal regimes after 2 months of incubation
# - metanalytical approach is used to quantify difference between predictions and 
#   measurements
# - the trajectory of freeze-thaw effect is extrapolated to later incubations
#===============================================================================
# ---- The same pipeline again ----
RootFT_obs <- TrLongerFull %>%
  filter(TempFactor %in% c("-5/5", "-18/5") & Origin == 'Roots') %>%
  mutate(Treatment = factor(TempFactor, levels = c("-5/5","-18/5")),
         Time = MonthsInc*30,
         Time_since_thaw = round((Time - 60.88) / 30.44),
         Locality = factor(Locality))

RootFT_obs %>%
  group_by(Treatment, Time_since_thaw) %>%
  summarize(mean_resp = mean(RootResp, na.rm=TRUE),
            se_resp = sd(RootResp, na.rm=TRUE)/sqrt(n()), .groups="drop")
# ---- Predictions ----
k_Pb_fun <- function(draws, temp) {
  exp(draws[["logk.(Intercept)"]] + draws[["logk.Temp"]] * temp)
}

a_Pb_fun <- function(draws, temp) {
  plogis(draws[["logit_a.(Intercept)"]] + draws[["logit_a.Temp"]] * temp)
}

predict_month2_root <- function(locality) {
  c0_val <- mean(RootFT_obs$C0[RootFT_obs$Locality == locality], na.rm=TRUE)
  
  k_neg5 <- k_Pb_fun(boot_coefs_Pb, -5)
  k_5    <- k_Pb_fun(boot_coefs_Pb, 5)
  a_neg5 <- a_Pb_fun(boot_coefs_Pb, -5)
  a_5    <- a_Pb_fun(boot_coefs_Pb, 5)
  
  pred_stable_neg5 <- k_neg5 * a_neg5 * c0_val * exp(-k_neg5 * 60.88)
  pred_stable_5    <- k_5    * a_5    * c0_val * exp(-k_5 * 60.88)
  
  tibble(
    Series = c("stable_-5","stable_5"),
    pred_mean = c(mean(pred_stable_neg5), mean(pred_stable_5)),
    pred_lwr  = c(quantile(pred_stable_neg5,0.025), quantile(pred_stable_5,0.025)),
    pred_upr  = c(quantile(pred_stable_neg5,0.975), quantile(pred_stable_5,0.975)),
    Locality = locality
  )
}

pred_val_root <- bind_rows(lapply(localities, predict_month2_root))
print(pred_val_root, n=20)

# ---- Observations ----
obs_month2_root <- TrRoots_cc %>%
  filter(Temp %in% c(-5,5), round((Time-60.88)/30.44)==0) %>%
  mutate(Series = ifelse(Temp==-5, "stable_-5","stable_5")) %>%
  group_by(Series, Locality) %>%
  summarize(obs_mean = mean(RootResp, na.rm=TRUE), n=n(), .groups="drop")

print(obs_month2_root, n=20)

# ---- formal test of reference point ----
compare_val_root <- pred_val_root %>%
  left_join(obs_month2_root, by = c("Series","Locality")) %>%
  mutate(within_band = obs_mean >= pred_lwr & obs_mean <= pred_upr)

print(compare_val_root, n=20)
compare_val_root %>% group_by(Series) %>% summarize(n_within = sum(within_band), n_total = n())

# ---- The validation was not entirely successful because unacknowledged locality effect in previous modelling ----
# ---- Re-validating across locality adjusted model calls ----
# Bootstraping both candidates
set.seed(123)
Pc_a_loc_boot <- boot_gnls(Pc_a_loc, nboot=999)
boot_coefs_Pc_a <- as.data.frame(Pc_a_loc_boot$t); names(boot_coefs_Pc_a) <- names(Pc_a_loc_boot$t0)
boot_coefs_Pc_a <- boot_coefs_Pc_a[complete.cases(boot_coefs_Pc_a),]

Pc_k_loc_boot <- boot_gnls(Pc_k_loc, nboot=999)
boot_coefs_Pc_k <- as.data.frame(Pc_k_loc_boot$t); names(boot_coefs_Pc_k) <- names(Pc_k_loc_boot$t0)
boot_coefs_Pc_k <- boot_coefs_Pc_k[complete.cases(boot_coefs_Pc_k),]

# Predictions
k_Pca_fun <- function(draws, temp) exp(draws[["logk.(Intercept)"]] + draws[["logk.Temp"]]*temp)
a_Pca_fun <- function(draws, locality, temp) {
  off <- if (locality==localities[1]) 0 else draws[[paste0("logit_a.Locality",locality)]]
  plogis(draws[["logit_a.(Intercept)"]] + off + draws[["logit_a.Temp"]]*temp)
}

k_Pck_fun <- function(draws, locality, temp) {
  off <- if (locality==localities[1]) 0 else draws[[paste0("logk.Locality",locality)]]
  exp(draws[["logk.(Intercept)"]] + off + draws[["logk.Temp"]]*temp)
}
a_Pck_fun <- function(draws, temp) plogis(draws[["logit_a.(Intercept)"]] + draws[["logit_a.Temp"]]*temp)

# re-validation against both
validate_model <- function(k_fun, a_fun, draws, has_a_loc) {
  bind_rows(lapply(localities, function(loc) {
    c0_val <- mean(TrRoots_cc$C0[TrRoots_cc$Locality == loc], na.rm=TRUE)
    k_neg5 <- if (has_a_loc) k_fun(draws, -5) else k_fun(draws, loc, -5)
    a_neg5 <- if (has_a_loc) a_fun(draws, loc, -5) else a_fun(draws, -5)
    pred <- k_neg5 * a_neg5 * c0_val * exp(-k_neg5*60.88)
    tibble(Locality=loc, pred_mean=mean(pred), pred_lwr=quantile(pred,0.025), pred_upr=quantile(pred,0.975))
  }))
}

val_Pca <- validate_model(k_Pca_fun, a_Pca_fun, boot_coefs_Pc_a, has_a_loc=TRUE) %>% left_join(obs_month2_root %>% filter(Series=="stable_-5"), by="Locality")
val_Pck <- validate_model(k_Pck_fun, a_Pck_fun, boot_coefs_Pc_k, has_a_loc=FALSE) %>% left_join(obs_month2_root %>% filter(Series=="stable_-5"), by="Locality")

val_Pca %>% mutate(within = obs_mean>=pred_lwr & obs_mean<=pred_upr)
val_Pck %>% mutate(within = obs_mean>=pred_lwr & obs_mean<=pred_upr)

# ---- re-validating across all ----
val_Pck_5 <- validate_model(k_Pck_fun, a_Pck_fun, boot_coefs_Pc_k, has_a_loc=FALSE) 
val_Pck_5 <- bind_rows(lapply(localities, function(loc) {
  c0_val <- mean(TrRoots_cc$C0[TrRoots_cc$Locality == loc], na.rm=TRUE)
  k5 <- k_Pck_fun(boot_coefs_Pc_k, loc, 5)
  a5 <- a_Pck_fun(boot_coefs_Pc_k, 5)
  pred <- k5 * a5 * c0_val * exp(-k5 * 60.88)
  tibble(Locality=loc, pred_mean=mean(pred), pred_lwr=quantile(pred,0.025), pred_upr=quantile(pred,0.975))
})) %>%
  left_join(obs_month2_root %>% filter(Series=="stable_5"), by="Locality") %>%
  mutate(within = obs_mean >= pred_lwr & obs_mean <= pred_upr)

print(val_Pck_5)

# ---- Quantifying the effect of freeze-thaw treatment ----
predict_5_5_root <- function(locality) {
  c0_val <- mean(TrRoots_cc$C0[TrRoots_cc$Locality == locality], na.rm=TRUE)
  k_neg5 <- k_Pck_fun(boot_coefs_Pc_k, locality, -5)
  k_5    <- k_Pck_fun(boot_coefs_Pc_k, locality, 5)
  a_5    <- a_Pck_fun(boot_coefs_Pc_k, 5)
  
  pred <- k_5 * a_5 * c0_val * exp(-(k_neg5 * 60.88 + k_5 * 14))
  tibble(pred_mean=mean(pred), pred_lwr=quantile(pred,0.025), pred_upr=quantile(pred,0.975))
}

pred_5_5_root_val <- bind_rows(lapply(localities, function(l) bind_cols(Locality=l, predict_5_5_root(l))))
print(pred_5_5_root_val)

predict_18_5_root <- function(locality) {
  c0_val <- mean(TrRoots_cc$C0[TrRoots_cc$Locality == locality], na.rm=TRUE)
  k_5 <- k_Pck_fun(boot_coefs_Pc_k, locality, 5)
  a_5 <- a_Pck_fun(boot_coefs_Pc_k, 5)
  
  pred <- k_5 * a_5 * c0_val * exp(-k_5 * 14)
  tibble(pred_mean=mean(pred), pred_lwr=quantile(pred,0.025), pred_upr=quantile(pred,0.975))
}

pred_18_5_root_val <- bind_rows(lapply(localities, function(l) bind_cols(Locality=l, predict_18_5_root(l))))
print(pred_18_5_root_val)

# ---- observations ----
obs_FT_root <- RootFT_obs %>%
  filter(TempFactor %in% c("-5/5","-18/5"), Time_since_thaw == 0) %>%
  group_by(Locality, TempFactor) %>%
  summarize(obs_mean = mean(RootResp, na.rm=TRUE), .groups="drop") %>%
  rename(Treatment = TempFactor)

print(obs_FT_root)

# ---- accounting for variance in measurements and predictions ----
# Observations
obs_summary_roots_5 <- RootFT_obs %>% filter(Treatment == "-5/5", Time_since_thaw == 0) %>%
  group_by(Locality) %>%
  summarize(obs_mean = mean(RootResp, na.rm=TRUE),
            obs_sd = sd(RootResp, na.rm=TRUE),
            n = n(), .groups = "drop")

# Calculating an excess with the proper variance 
meta_data_roots_5 <- pred_5_5_root_val %>%   # control-soil predicted values, from the earlier N2-based check
  left_join(obs_summary_roots_5, by = "Locality") %>%
  mutate(
    obs_se = obs_sd / sqrt(n),
    pred_se = (pred_upr - pred_lwr) / (2*1.96),
    yi = log(obs_mean) - log(pred_mean),
    vi = (obs_se/obs_mean)^2 + (pred_se/pred_mean)^2
  )

print(meta_data_roots_5)

# Formal test
res_overall_roots_5 <- rma(yi, vi, data = meta_data_roots_5, method = "REML")
summary(res_overall_roots_5)

# back-transformed pooled excess ratio
exp(c(estimate = coef(res_overall_roots_5), ci.lb = res_overall_roots_5$ci.lb, ci.ub = res_overall_roots_5$ci.ub))

# try elevation modifier
meta_data_roots_5_elev <- meta_data_roots_5 %>% left_join(elev_lookup, by="Locality")
res_roots_5_elev <- rma(yi, vi, mods = ~ Elevation, data = meta_data_roots_5_elev, method="REML")
summary(res_roots_5_elev)
res_roots_5_veg <- rma(yi, vi, mods = ~ Vegetation, data = meta_data_roots_5_elev, method="REML")
summary(res_roots_5_veg)

# Observations
obs_summary_roots_18 <- RootFT_obs %>% filter(Treatment == "-18/5", Time_since_thaw == 0) %>%
  group_by(Locality) %>%
  summarize(obs_mean = mean(RootResp, na.rm=TRUE),
            obs_sd = sd(RootResp, na.rm=TRUE),
            n = n(), .groups = "drop")

# Calculating an excess with the proper variance 
meta_data_roots_18 <- pred_18_5_root_val %>%   # control-soil predicted values, from the earlier N2-based check
  left_join(obs_summary_roots_18, by = "Locality") %>%
  mutate(
    obs_se = obs_sd / sqrt(n),
    pred_se = (pred_upr - pred_lwr) / (2*1.96),
    yi = log(obs_mean) - log(pred_mean),
    vi = (obs_se/obs_mean)^2 + (pred_se/pred_mean)^2
  )

print(meta_data_roots_18)

# Formal test
res_overall_roots_18 <- rma(yi, vi, data = meta_data_roots_18, method = "REML")
summary(res_overall_roots_18)

# back-transformed pooled excess ratio
exp(c(estimate = coef(res_overall_roots_18), ci.lb = res_overall_roots_18$ci.lb, ci.ub = res_overall_roots_18$ci.ub))

# try elevation modifier
meta_data_roots_18_elev <- meta_data_roots_18 %>% left_join(elev_lookup, by="Locality")
res_roots_18_elev <- rma(yi, vi, mods = ~ Elevation, data = meta_data_roots_18_elev, method="REML")
summary(res_roots_18_elev)
res_roots_18_veg <- rma(yi, vi, mods = ~ Vegetation, data = meta_data_roots_18_elev, method="REML")
summary(res_roots_18_veg)

#==============time to excess respiration attenuation
set.seed(123)
n_boot <- 2000

# 1. Bootstrap draws of the month-0 excess ratio, from the rma() model's own sampling distribution
ratio0_draws_roots <- exp(rnorm(n_boot, mean = coef(res_overall_roots_5), sd = res_overall_roots_5$se))

# 2. Bootstrap the relative-decline trajectory by resampling replicates within each Series/Time cell
rel_boot_roots <- replicate(n_boot, {
  boot_sample <- rel_decline_root_derived_clean %>% filter(TempFactor == "-5/5") %>%
    group_by(Locality, Time_since_thaw) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    summarize(mean_rel = mean(rel_change, na.rm=TRUE), .groups="drop")
  boot_sample$mean_rel
}) %>% t()   # n_boot x 3 matrix, columns = Time_since_thaw 0,1,2

colnames(rel_boot_roots) <- unlist(lapply(unique(TrLongerFull$Locality), paste0, c("-m0", "-m2", "-m3")))


# 3. Combine: ratio(m) = ratio(0) * rel_obs(m)/rel_obs(0), but rel_obs(0) = 1 by construction
ratio_traj_boot_roots <- ratio0_draws_roots * rel_boot_roots   # element-wise: each row is one bootstrap replicate's full trajectory

ratio_traj_summary_roots <- tibble(
  Legend = colnames(rel_boot_roots),
  median = apply(ratio_traj_boot_roots, 2, median, na.rm = T),
  lwr = apply(ratio_traj_boot_roots, 2, quantile, 0.025, na.rm = T),
  upr = apply(ratio_traj_boot_roots, 2, quantile, 0.975, na.rm = T)
)
print(ratio_traj_summary)

# boot_coefs_Pc_k and ratio0_draws/rel_boot have different lengths -> resample to align
# (treats the two models' bootstraps as independent)
set.seed(123)

k_Pck_fun <- function(draws, locality, temp) {
  off <- if (locality==localities[1]) 0 else draws[[paste0("logk.Locality",locality)]]
  exp(draws[["logk.(Intercept)"]] + off + draws[["logk.Temp"]]*temp)
}

idx <- sample(nrow(boot_coefs_Pc_k), n_boot, replace = TRUE)

kn_roots <- sapply(unique(TrLongerFull$Locality), k_Pck_fun, simplify = 'matrix', draws = boot_coefs_Pc_k[idx, ], temp = -5)
kn_roots <- kn_roots[, rep(1:6, each = 3)]
kn_roots[, seq(1, 18, by = 3)] <- c(1)
kn_roots[, seq(2, 18, by = 3)] <- exp(-kn_roots[, seq(2, 18, by = 3)]*30)
kn_roots[, seq(3, 18, by = 3)] <- exp(-kn_roots[, seq(3, 18, by = 3)]*30*2)

kp_roots <- sapply(unique(TrLongerFull$Locality), k_Pck_fun, simplify = 'matrix', draws = boot_coefs_Pc_k[idx, ], temp = 5)
kp_roots <- kp_roots[, rep(1:6, each = 3)]
kp_roots[, seq(1, 18, by = 3)] <- c(1)
kp_roots[, seq(2, 18, by = 3)] <- exp(-kp_roots[, seq(2, 18, by = 3)]*30)
kp_roots[, seq(3, 18, by = 3)] <- exp(-kp_roots[, seq(3, 18, by = 3)]*30*2)


scen_roots <- list(all_cold = kn_roots,
                   all_warm = kp_roots)

traj_roots <- lapply(scen_roots, function(e) ratio0_draws_roots * rel_boot_roots / e)

trajDFrame_roots <- data.frame(Scenario = rep(c("Cold", "Warm"), each = 2000*18),
                               Locality = rep(rep(unique(TrLongerFull$Locality), each = 3*2000), 2),
                               Month = rep(rep(c(0, 1, 2), each = 2000), times = 6*2),
                               Rm = as.numeric(c(traj_roots$all_cold, traj_roots$all_warm))) 

trajDFrame_roots |> group_by(Scenario, Month) |> filter(Rm < 1000) |> 
  summarise(Center = mean(Rm, na.rm = T),
            lb = quantile(Rm, probs = c(0.025), na.rm = T),
            ub = quantile(Rm, probs = c(0.975), na.rm = T))

#==============time to excess respiration attenuation
set.seed(123)
n_boot <- 2000

# 1. Bootstrap draws of the month-0 excess ratio, from the rma() model's own sampling distribution
ratio0_draws_roots18 <- exp(rnorm(n_boot, mean = coef(res_overall_roots_18), sd = res_overall_roots_18$se))

# 2. Bootstrap the relative-decline trajectory by resampling replicates within each Series/Time cell
rel_boot_roots18 <- replicate(n_boot, {
  boot_sample <- rel_decline_root_derived_clean %>% filter(TempFactor == "-18/5") %>%
    group_by(Locality, Time_since_thaw) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    summarize(mean_rel = mean(rel_change, na.rm=TRUE), .groups="drop")
  boot_sample$mean_rel
}) %>% t()   # n_boot x 3 matrix, columns = Time_since_thaw 0,1,2

colnames(rel_boot_roots18) <- unlist(lapply(unique(TrLongerFull$Locality), paste0, c("-m0", "-m2", "-m3")))


# 3. Combine: ratio(m) = ratio(0) * rel_obs(m)/rel_obs(0), but rel_obs(0) = 1 by construction
ratio_traj_boot_roots18 <- ratio0_draws_roots18 * rel_boot_roots18   # element-wise: each row is one bootstrap replicate's full trajectory

ratio_traj_summary_roots18 <- tibble(
  Legend = colnames(rel_boot_roots18),
  median = apply(ratio_traj_boot_roots18, 2, median, na.rm = T),
  lwr = apply(ratio_traj_boot_roots18, 2, quantile, 0.025, na.rm = T),
  upr = apply(ratio_traj_boot_roots18, 2, quantile, 0.975, na.rm = T)
)
print(ratio_traj_summary)

# boot_coefs_N2 and ratio0_draws/rel_boot have different lengths -> resample to align
scen_roots18 <- list(all_cold = matrix(rep(1, times = 2000*18), ncol = 18),
                     all_warm = kp_roots)

traj_roots18 <- lapply(scen_roots18, function(e) ratio0_draws_roots18 * rel_boot_roots18 / e)

trajDFrame_roots18 <- data.frame(Scenario = rep(c("Cold", "Warm"), each = 2000*18),
                                 Locality = rep(rep(unique(TrLongerFull$Locality), each = 3*2000), 2),
                                 Month = rep(rep(c(0, 1, 2), each = 2000), times = 6*2),
                                 Rm = as.numeric(c(traj_roots18$all_cold, traj_roots18$all_warm))) 

trajDFrame_roots18 |> group_by(Scenario, Month) |> summarise(Center = mean(Rm, na.rm = T),
                                                             lb = quantile(Rm, probs = c(0.025), na.rm = T),
                                                             ub = quantile(Rm, probs = c(0.975), na.rm = T))

# ---- Direct comparison of measurements at month 4 ----
#1. -5/5
anova(aov(RootResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Origin == 'Roots')))
anova(aov(RootResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5' & Origin == 'Roots')))
TukeyHSD(aov(RootResp ~ TempFactor*Locality, data = TrLongerFull,
             subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5'  & Roots == 'Yes')))

Root5end <- as.data.frame(TukeyHSD(aov(RootResp ~ TempFactor*Locality, data = TrLongerFull,
                                       subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-18/5'  & Roots == 'Yes')))$`TempFactor:Locality`)
Root5end$comp <- rownames(Root5end)

for(i in unique(TrLonger$Locality)){
  #Soil specific comparison
  comp5negative <- paste0('-5/5:', i, '--5:', i)
  comp5positive <- paste0('5:', i, '--5/5:', i)
  
  print(rbind(
    Root5end[grepl(comp5negative, Root5end$comp), c("diff", 'p adj', 'comp')],
    Root5end[grepl(comp5positive, Root5end$comp), c("diff", 'p adj', 'comp')]
  ))
  
}

#1. -18/5
anova(aov(SoilResp ~ TempFactor, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5'  & Origin == 'Roots')))
anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull,
          subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5'  & Origin == 'Roots')))
TukeyHSD(aov(SoilResp ~ TempFactor, data = TrLongerFull,
             subset = (MonthsInc == 4 & TempFactor != '0' & TempFactor != '-5/5'  & Origin == 'Roots')))



#===============================================================================
# Figures for root-amended soils - SOC-derived respiration
#===============================================================================
#=======================Fraction of root-derived respiration in respect to total respiration for constant thermal regimes
RootsFraction_StableFig <- TrRoots_cc %>%
  mutate(Locality = factor(Locality,
                           levels = c("SC", "VW", "NW", "FU", "JA", "RA"),
                           labels = c(expression(atop("2160 m", "scree")),
                                      expression(atop("2160 m", "meadow")),
                                      expression(atop("2130 m", "meadow")),
                                      expression(atop("1780 m", "meadow")),
                                      expression(atop("1470 m", "forest")),
                                      expression(atop("1330 m", "forest"))))) |> 
  group_by(Locality, Temp, MonthsInc) |> summarise(y = mean(frRoot, na.rm = T),
                                                   ySD = sd(frRoot, na.rm = T)) |> 
  ggplot(aes(MonthsInc, y)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 1), aes(fill = factor(Temp)),
           color = 'black') + facet_wrap(. ~ Locality,
                                         nrow = 2, labeller = label_parsed) +
  mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD, color = factor(Temp)), 
                          position = position_dodge(width = 1), show.legend = F) +
  scale_fill_manual(values = c('white', 'grey', 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  labs(color = 'Thermal regime', fill = "Thermal regime", x = 'Time (months)',
       y = expression(paste("Root-derived respiration (nondimensional fraction)"))) +
  theme(legend.position = c(0.08, 0.36))

ggsave("StatOut/RootsFraction_StableFig.pdf", RootsFraction_StableFig,
       width = 300, height = 220, units = "mm", device = cairo_pdf)

#=================Measured and predicted respiration rates
#Measured
RootStable <- TrRoots_cc |> group_by(Locality, Temp, Time) |> summarise(y = mean(RootResp, na.rm = T),
                                                                        ySD = sd(RootResp, na.rm = T),
                                                                        TOC_anchor = median(C0, na.rm = T))
#Predicted
RootStable$yhat <- numeric(length = nrow(RootStable))
RootStable$lwr <- numeric(length = nrow(RootStable))
RootStable$upr <- numeric(length = nrow(RootStable))

for(i in 1:nrow(RootStable)){
  #Model parameters for row i
  loc_i <- RootStable$Locality[i]
  if(loc_i == 'FU'){
    k_i = exp(boot_coefs_Pc_k[['logk.(Intercept)']] + boot_coefs_Pc_k[['logk.Temp']]*RootStable$Temp[i])
    
  } else {
    k_i = exp(boot_coefs_Pc_k[['logk.(Intercept)']] + 
                boot_coefs_Pc_k[[paste0('logk.Locality', loc_i)]] +
                boot_coefs_Pc_k[['logk.Temp']]*RootStable$Temp[i])
  }
  
  alfa_i =  plogis(boot_coefs_Pc_k[['logit_a.(Intercept)']] + 
                     boot_coefs_Pc_k[['logit_a.Temp']]*RootStable$Temp[i])
  
  
  
  #Prediction for row i
  RootStable[i, c('yhat', 'lwr', 'upr')] <- matrix(c(mean(RootStable$TOC_anchor[i]*
                                                            k_i*alfa_i*exp(-k_i*RootStable$Time[i])), 
                                                     quantile(x = RootStable$TOC_anchor[i]*
                                                                k_i*alfa_i*exp(-k_i*RootStable$Time[i]),
                                                              probs = c(0.025, 0.975))), nrow = 1)
  
  
}

RootAmendedStable_RootderivedFig <- RootStable %>%
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
       y = expression(paste("Root-derived respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
  new_scale_fill() +
  new_scale_color() +
  geom_point(cex = 4, pch = 21, aes(y = yhat, fill = factor(Temp)), position = position_dodge(width = 1),
             show.legend = F) + 
  geom_errorbar(aes(ymin = lwr, ymax = upr, colour = factor(Temp)), 
                position = position_dodge(width = 1), width = 0.1, show.legend = F) +
  scale_fill_manual(values = c('red', 'red', 'red')) +
  scale_color_manual(values = c('red', 'red', 'red')) +
  theme(legend.position = c(0.08, 0.88))

ggsave("StatOut/RootAmendedStable_RootderivedFig.pdf", RootAmendedStable_RootderivedFig,
       width = 300, height = 180, units = "mm", device = cairo_pdf)

#=================Predicted respiration rates - continuous time dependent change in constant thermal regimes
#Simulations are stored here
RootContinuous <- tibble(Locality = rep(unique(RootStable$Locality), each = 25*3),
                         Time = rep(seq(0, 120, by = 5), times = 6*3),
                         Temp = rep(rep(c(-5, 0, 5), each = 25), times = 6),
                         yhat = numeric(length = length(Locality)),
                         lwr = numeric(length = length(Locality)),
                         upr = numeric(length = length(Locality)))
RootContinuous <- TrRoots_cc |> group_by(Locality) |> summarise(C0 = mean(C0), .groups = "drop") |> 
  merge(RootContinuous, by = "Locality")

for(i in 1:nrow(RootContinuous)){
  #Model parameters for row i
  loc_i <- RootContinuous$Locality[i]
  if(loc_i == 'FU'){
    k_i = exp(boot_coefs_Pc_k[['logk.(Intercept)']] + boot_coefs_Pc_k[['logk.Temp']]*RootContinuous$Temp[i])
    
  } else {
    k_i = exp(boot_coefs_Pc_k[['logk.(Intercept)']] + 
                boot_coefs_Pc_k[[paste0('logk.Locality', loc_i)]] +
                boot_coefs_Pc_k[['logk.Temp']]*RootContinuous$Temp[i])
  }
  
  alfa_i =  plogis(boot_coefs_Pc_k[['logit_a.(Intercept)']] + 
                     boot_coefs_Pc_k[['logit_a.Temp']]*RootContinuous$Temp[i])
  
  
  
  #Prediction for row i
  RootContinuous[i, c('yhat', 'lwr', 'upr')] <- matrix(c(mean(RootContinuous$C0[i]*
                                                            k_i*alfa_i*exp(-k_i*RootContinuous$Time[i])), 
                                                     quantile(x = RootContinuous$C0[i]*
                                                                k_i*alfa_i*exp(-k_i*RootContinuous$Time[i]),
                                                              probs = c(0.05, 0.95))), nrow = 1)
  
  
}

RootSimulation_Fig <- 
  RootContinuous %>%
  mutate(Locality = factor(Locality,
                           levels = c("SC", "VW", "NW", "FU", "JA", "RA"),
                           labels = c(expression(atop("2160 m", "scree")),
                                      expression(atop("2160 m", "meadow")),
                                      expression(atop("2130 m", "meadow")),
                                      expression(atop("1780 m", "meadow")),
                                      expression(atop("1470 m", "forest")),
                                      expression(atop("1330 m", "forest"))))) |> 
  ggplot(aes(Time/30, yhat)) + 
  geom_point(pch = 21, cex = 4, aes(fill = factor(Temp))) + 
  facet_wrap(. ~ Locality, nrow = 2, labeller = label_parsed) +
  mytheme + geom_line(aes(color = factor(Temp))) +
  scale_fill_manual(values = c('white', 'grey', 'black')) +
  scale_color_manual(values = c('black', 'black', 'black')) +
  labs(color = 'Thermal regime', fill = "Thermal regime", x = 'Time (months)',
       y = expression(paste("Root-derived respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
  theme(legend.position = c(0.08, 0.68)) + scale_y_log10(limits = c(0.005, 100), breaks = c(0.01, 0.1, 1, 10, 100),
                                                        labels = c(0.01, 0.1, 1, 10, 100)) +
    geom_vline(xintercept = 2, col = 'blue') + geom_vline(xintercept = 4, col = 'red')

ggsave("StatOut/RootSimulation_Fig.pdf", RootSimulation_Fig,
       width = 320, height = 200, units = "mm", device = cairo_pdf)

#=======================Fraction of root-derived respiration in respect to total respiration for alternating thermal regimes
RootsFraction_AlternatingFig <- TrLongerFull %>% filter(Origin == 'Roots' & StabTemp == 'No') |> 
  mutate(Locality = factor(Locality,
                           levels = c("SC", "VW", "NW", "FU", "JA", "RA"),
                           labels = c(expression(atop("2160 m", "scree")),
                                      expression(atop("2160 m", "meadow")),
                                      expression(atop("2130 m", "meadow")),
                                      expression(atop("1780 m", "meadow")),
                                      expression(atop("1470 m", "forest")),
                                      expression(atop("1330 m", "forest")))),
         Temp = factor(Temp, levels = c('-5', '-18'),
                       labels = c("-5˚C ⟶ +5˚C", "-18˚C ⟶ +5˚C"))) |> 
  group_by(Locality, Temp, MonthsInc) |> summarise(y = mean(frRoot, na.rm = T),
                                                   ySD = sd(frRoot, na.rm = T)) |> 
  ggplot(aes(MonthsInc, y)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 1), aes(fill = Temp),
           color = 'black') + facet_wrap(. ~ Locality,
                                         nrow = 2, labeller = label_parsed) +
  mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD, color = Temp), 
                          position = position_dodge(width = 1), show.legend = F) +
  scale_fill_manual(values = c('dodgerblue1', 'blue4')) +
  scale_color_manual(values = c('black', 'black')) +
  labs(color = 'Thermal regime', fill = "Thermal regime", x = 'Time (months)',
       y = expression(paste("Root-derived respiration (nondimensional fraction)"))) +
  theme(legend.position = c(0.1, 0.36))

ggsave("StatOut/RootsFraction_AlternatingFig.pdf", RootsFraction_AlternatingFig,
       width = 300, height = 220, units = "mm", device = cairo_pdf)

#=======================Freeze-thaw cycles - composite figure
pred_month2_FT_root <- bind_rows(pred_5_5_root_val, pred_18_5_root_val) |> 
  mutate(TempFactor = rep(c("-5/5", "-18/5"), each = 6)) 


RootRef <- TrRoots_cc |> group_by(Locality, Temp, MonthsInc) |> 
  summarise(yref = mean(RootResp, na.rm = T), yrefSD = sd(RootResp, na.rm = T)) |> 
  filter(Temp != 0) |> filter(MonthsInc != 3) |>  
  mutate(TempF = factor(Temp, levels = c(-5, 5), labels = c("-5˚C", '+5˚C')))

pred_month2_FT_root2 <- rbind(pred_month2_FT_root, pred_month2_FT_root)
pred_month2_FT_root2$MonthsInc <- rep(c(2, 4), each = nrow(pred_month2_FT_root)) 
pred_month2_FT_root2[(nrow(pred_month2_FT_root)+1):nrow(pred_month2_FT_root2), c(2:4)] <- NA

RootDerivedFT_compositeFig <- TrLongerFull |> filter(Origin == 'Roots' & StabTemp == 'No' & MonthsInc != 3) |> 
  group_by(Locality, TempFactor, MonthsInc) |> summarise(y = mean(RootResp, na.rm = T),
                                                         ySD = sd(RootResp, na.rm = T)) |> 
  inner_join(pred_month2_FT_root2, by = c('Locality', "TempFactor", "MonthsInc")) |> 
  full_join(RootRef) |> 
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
       y = expression(paste("Root-derived respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
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

ggsave("StatOut/RootDerivedFT_compositeFig.pdf", RootDerivedFT_compositeFig,
       width = 340, height = 180, units = "mm", device = cairo_pdf)

#==================Trajectory of freeze-thaw effect on root-derived respiration rates in root-amended soils 
trajDFrame_roots$Series <- c('-5/5')
trajDFrame_roots18$Series <- c('-18/5')

RootDerivedFT_temporalFig <- rbind(trajDFrame_roots, trajDFrame_roots18) |> group_by(Series, Scenario, Locality, Month) |> 
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
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 16)) +
  labs(y = c("Relative effect of alternating thermal regime (nondimensional)"),
       color = c("Time (months)"), fill = c("Time (months)"))

ggsave("StatOut/RootDerivedFT_temporalFig.pdf", RootDerivedFT_temporalFig,
       width = 360, height = 180, units = "mm", device = cairo_pdf)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF SCRIPT