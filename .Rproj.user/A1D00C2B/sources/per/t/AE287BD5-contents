####################################################################################################
###                        Tatra Mountains soil cores incubation experiment                      ###
######                                      Statistical analysis                             #######
####################################################################################################
#===============================Libraries + ggplot theme
mytheme <- readRDS('/mnt/9C64E5F864E5D554/Documents and Settings/capekp00/Documents/helpfull_R_Matlab_script/ggtheme.rds')
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyr)
library(nlme)
library(minpack.lm)
library(nlraa)
library(car)
#==============================END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#===============================Data
TrLongerFull <- read.csv('TrLongerFull.csv')

#==============================END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
####################################################################################################
#========================================Statistical analysis======================================#
####################################################################################################
#===========================1. Stable temperatures (-5, 0, 5) + soil respiration only [non-linear mixed effects]
##=================1a. Control soils - nested models
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

# ---- correlation of C0 with total soil organic carbon - there should be thight correlation ----
# reconstruct absolute per-soil C0 (intercept + offset, 0 for reference level FU)
locality_levels <- levels(Ctr$Locality)
coefs <- coef(M1_gnls)

C0_abs <- c(coefs["C0.(Intercept)"],
            coefs["C0.(Intercept)"] + coefs[grep("^C0\\.Locality", names(coefs))])
names(C0_abs) <- locality_levels

C0_df <- tibble(Locality = locality_levels, C0 = C0_abs)

# mean TOC per soil, pooled across all temps (control only, as discussed)
toc_mean_df <- Ctr %>%
  group_by(Locality) %>%
  summarize(TOC_mean = max(soil_C, na.rm = TRUE))

compare_df <- left_join(C0_df, toc_mean_df, by = "Locality")

ggplot(compare_df, aes(TOC_mean, C0, label = Locality)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8) +
  labs(x = "Mean TOC (soil, all temps pooled)", y = "Fitted C0 (active pool)")

# ---- Is the Gaussian temperature sensitivity function more accurate representation? ----
M1_quad <- update(M1_gnls, 
                  params = list(C0 ~ Locality, logk ~ Temp + I(Temp^2)),
                  start = c(C0 = c(1824.16, -76.99, -1019.85, -791.49, -1376.17, -946.77),
                            logk = c(-5.7250, 0.1758, 0)))

anova(M1_gnls, M1_quad) #It is not !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#3. The decomposition constant k differs among soils - it is plausible but uncertain ----
M2_gnls <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = Ctr,
  params = list(C0 ~ Locality, logk ~ Locality + Temp),
  start = c(C0 = c(1824.16, -76.99, -1019.85, -791.49, -1376.17, -946.77),
            logk = c(-5.7250, 0, 0, 0, 0, 0, 0.1758)),
  control = gnlsControl(maxIter = 300, msMaxIter = 300),
  na.action = na.exclude
)

summary(M2_gnls)
anova(M1_gnls, M2_gnls)
AIC(M1_gnls, M2_gnls) #There is no significant difference !!!!!!!!!!!!!!!!!!

# ---- fitting heteroscedascity on M1 before moving forward ----
M1_diag <- tibble(fitted = fitted(M1_gnls), resid = residuals(M1_gnls),
                  Temp = Ctr$Temp, Locality = Ctr$Locality)

ggplot(M1_diag, aes(fitted, resid)) +
  geom_point(alpha=0.5) + geom_hline(yintercept=0, linetype="dashed") +
  labs(title = "Residuals vs fitted, M1")

ggplot(M1_diag, aes(factor(Temp), abs(resid))) + geom_boxplot() +
  labs(title = "Absolute residuals by Temp, M1")

M1_varIdent_temp <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = Ctr_cc,
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


#============================Refitting with known SOC and fraction of available TOC instead of C0
toc_anchor_soil <- Ctr %>% group_by(Locality) %>%
  summarize(TOC_anchor = max(soil_C*1000, na.rm=TRUE))   # pooled across all temp/FT, same units as SoilResp

Ctr_cc <- Ctr_cc %>% left_join(toc_anchor_soil, by = "Locality")

# a=1 baseline vs free shared a vs a~Temp
M1_a1   <- M1_varIdent_temp   # existing model, for AIC reference

M1_aB <- gnls(
  SoilResp ~ exp(logk) * plogis(logit_a) * TOC_anchor * exp(-exp(logk) * Time),
  data = Ctr_cc,
  params = list(logk ~ Temp, logit_a ~ 1),
  start = c(logk = coef(M1_gnls)[c("logk.(Intercept)","logk.Temp")], logit_a = qlogis(0.01)),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)

M1_aC <- update(M1_aB, params = list(logk ~ Temp, logit_a ~ Temp),
                start = c(coef(M1_aB), 0))

anova(M1_aB, M1_aC)   # does available fraction itself depend on temperature?

logit_a_int  <- coef(M1_aC)["logit_a.(Intercept)"]
logit_a_temp <- coef(M1_aC)["logit_a.Temp"]

a_at <- function(T) plogis(logit_a_int + logit_a_temp * T)

tibble(Temp = c(-5, 0, 5), a = sapply(c(-5,0,5), a_at))

summary(lm(MBC ~ Temp, Ctr_cc))

##=================1a. Control soils - nested models END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##=================1b. Control soils and soils with roots - nested models
TrSoil <- TrLongerFull |> filter(Origin != "Roots" & StabTemp == 'Yes')
TrSoil$Locality <- as.factor(TrSoil$Locality)
TrSoil$Origin <- as.factor(TrSoil$Origin)
# ---- Convert months of incubation to days of incubation to match SoilResp in ug/d/d ----
TrSoil$Time <- TrSoil$MonthsInc*30
nrow(TrSoil)          # should be 324 (minus any NAs, same as before)
table(TrSoil$Origin)   

# ---- nested non-linear mixed effect models - appending to M1_gnls ----
#4. Refitting M1_gnls with new dataset that include respiration of soil organic carbon from soils exposed to roots
R0 <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = TrSoil,
  params = list(C0 ~ Locality, logk ~ Temp),
  start = c(C0 = c(1824.16, -76.99, -1019.85, -791.49, -1376.17, -946.77),
            logk = c(-5.7250, 0.1758)),
  control = gnlsControl(maxIter = 300, msMaxIter = 300),
  na.action = na.exclude
)
summary(R0)
cbind(coef(M1_gnls), coef(R0))

#5. Does the root addition speeds the decay rate up, i.e. increases k uniformly across localities and temperatures?
R1 <- update(R0, params = list(C0 ~ Locality, logk ~ Temp + Origin),
             start = c(coef(R0)[1:6], coef(R0)[7:8], 0))
#6. Does the root addition changes the temperature sensitivity of soil organic carbon decay rate?
R2 <- update(R1, params = list(C0 ~ Locality, logk ~ Temp * Origin),
             start = c(coef(R1), 0))
#7. Does the root addition also change pool size C0 on top of any kinetic effect?
R3 <- update(R2, params = list(C0 ~ Locality + Origin, logk ~ Temp * Origin),
             start = c(coef(R2)[1:6], 0, coef(R2)[7:10]))

anova(R0, R1, R2, R3)
AIC(R0, R1, R2, R3)

# ---- the best model is R2. The specific information will be extracted from it ----
summary(R2)
anova(R2)
# ---- Bootstraping R2 to calculate confidence intervals of parameters and predictions ----
set.seed(123)
# model must be refit to dataframe without NAs to make bootstraping working
TrSoil_cc <- TrSoil[complete.cases(TrSoil[, c("SoilResp", "Temp", "Time", "Locality", "Origin")]), ]
nrow(TrSoil_cc)   # should be 288

R2_refit <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = TrSoil_cc,
  params = list(C0 ~ Locality, logk ~ Temp * Origin),
  start = coef(R2),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
  # no na.action needed — there are no NAs left to handle
)

length(fitted(R2_refit))   # should now be exactly 288, no padding
identical(round(coef(R2_refit), 4), round(coef(R2), 4))   # confirm same fit

R2_boot <- boot_gnls(R2_refit, nboot = 999)
str(R2_boot)

# ---- modelling heteroscedascity
R2_varIdent_temp <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = TrSoil_cc,
  params = list(C0 ~ Locality, logk ~ Temp * Origin),
  start = coef(R2_refit),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)

anova(R2_refit, R2_varIdent_temp)
summary(R2_varIdent_temp)
cbind(R2_refit = coef(R2_refit), R2_varIdent_temp = coef(R2_varIdent_temp))

# ---- Plotting predictions against observations ----
TrSoil_fit_check <- TrSoil_cc %>%   # the actual data R2 was fit on
  mutate(predicted = fitted(R2_varIdent_temp),   # or predict(R2_refit) if fitted() doesn't apply cleanly to gnls
         observed = SoilResp)

ggplot(TrSoil_fit_check, aes(predicted, observed, color = interaction(Temp, Origin))) +
  geom_abline(slope = 1, intercept = 0, color = "red", lwd = 1.5) +
  geom_point(size = 2, alpha = 0.6) +
  mytheme +
  labs(x = "Predicted respiration", y = "Observed respiration")

ggplot(TrSoil_fit_check, aes(predicted, observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  geom_point(aes(color = factor(Time)), size = 2, alpha = 0.6) +
  coord_fixed() +
  labs(color = "Time (days)", title = "R2 fit quality by measurement time")
# ---- rebooting the improved model ----
set.seed(123)
R2_varIdent_boot <- boot_gnls(R2_varIdent_temp, nboot = 999)
boot_coefs_v2 <- as.data.frame(R2_varIdent_boot$t)
names(boot_coefs_v2) <- names(R2_varIdent_boot$t0)
boot_coefs_v2 <- boot_coefs_v2[complete.cases(boot_coefs_v2), ]


# ---- extracting coefficients with uncertainty estimates ----
b0  <- boot_coefs_v2[["logk.(Intercept)"]]
bT  <- boot_coefs_v2[["logk.Temp"]]
bR  <- boot_coefs_v2[["logk.OriginSoil with roots"]]
bTR <- boot_coefs_v2[["logk.Temp:OriginSoil with roots"]]

# ---- k-vs-temperature curve for soils with and without the root addition, with 95% bootstrap intervals ----
temp_grid <- seq(-5, 5, by = 0.5)

k_curve_df <- lapply(temp_grid, function(T) {
  logk_absent  <- b0 + bT * T
  logk_present <- b0 + bT * T + bR + bTR * T
  tibble(
    Temp = T,
    Root = c("absent", "present"),
    k_mean = c(mean(exp(logk_absent)), mean(exp(logk_present))),
    k_lwr  = c(quantile(exp(logk_absent), 0.025), quantile(exp(logk_present), 0.025)),
    k_upr  = c(quantile(exp(logk_absent), 0.975), quantile(exp(logk_present), 0.975))
  )
}) %>% bind_rows()

print(k_curve_df)

# ---- Q10 comparison, with a properly-correlated difference interval ----
Q10_absent  <- exp(10 * bT)
Q10_present <- exp(10 * (bT + bTR))
Q10_diff    <- Q10_present - Q10_absent

Q10_df <- tibble(
  Root = c("absent", "present", "difference"),
  Q10_mean = c(mean(Q10_absent), mean(Q10_present), mean(Q10_diff)),
  Q10_lwr  = c(quantile(Q10_absent, 0.025), quantile(Q10_present, 0.025), quantile(Q10_diff, 0.025)),
  Q10_upr  = c(quantile(Q10_absent, 0.975), quantile(Q10_present, 0.975), quantile(Q10_diff, 0.975))
)

print(Q10_df)

# ---- Rescale to common temperature gradient ----
R_gas <- 8.314  # J/mol/K

T1_cold <- 268.15; T2_cold <- 278.15   # -5 to 5 °C
T1_warm <- 298.15; T2_warm <- 308.15   # 25 to 35 °C

bracket_cold <- 1/T1_cold - 1/T2_cold
bracket_warm <- 1/T1_warm - 1/T2_warm

# Ea per bootstrap draw, for each root treatment (reusing Q10_absent/Q10_present from before)
Ea_absent  <- R_gas * log(Q10_absent)  / bracket_cold
Ea_present <- R_gas * log(Q10_present) / bracket_cold

# recalculated Q10 at 25-35°C, same Ea, same bootstrap draws
Q10_warm_absent  <- exp(Ea_absent  / R_gas * bracket_warm)
Q10_warm_present <- exp(Ea_present / R_gas * bracket_warm)

Q10_rescaled_df <- tibble(
  Root = c("absent", "present"),
  Ea_kJmol = c(mean(Ea_absent), mean(Ea_present)) / 1000,
  Q10_cold_mean = c(mean(Q10_absent), mean(Q10_present)),
  Q10_warm_mean = c(mean(Q10_warm_absent), mean(Q10_warm_present)),
  Q10_warm_lwr  = c(quantile(Q10_warm_absent, 0.025), quantile(Q10_warm_present, 0.025)),
  Q10_warm_upr  = c(quantile(Q10_warm_absent, 0.975), quantile(Q10_warm_present, 0.975))
)

print(Q10_rescaled_df)

# ---- plotting the R2 results ----
ggplot(k_curve_df, aes(Temp, k_mean, color = Root, fill = Root)) +
  geom_ribbon(aes(ymin = k_lwr, ymax = k_upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Temperature (°C)", y = "Kinetic constant k (1/d)") +
  mytheme + theme(legend.position = c(0.15, 0.85))

ggplot(Q10_df %>% filter(Root != "difference"), aes(Root, Q10_mean)) +
  geom_pointrange(aes(ymin = Q10_lwr, ymax = Q10_upr), linewidth = 1) +
  labs(y = expression(Q[10]), x = NULL) + mytheme

# ---- Calculating the crossover temperature at which priming sign flips ----
bR  <- boot_coefs_v2[["logk.OriginSoil with roots"]]
bTR <- boot_coefs_v2[["logk.Temp:OriginSoil with roots"]]

T_cross <- -bR / bTR

tibble(
  metric = "T_cross (priming sign flip)",
  mean = mean(T_cross, na.rm=TRUE),
  median = median(T_cross, na.rm=TRUE),
  lwr = quantile(T_cross, 0.025, na.rm=TRUE),
  upr = quantile(T_cross, 0.975, na.rm=TRUE)
)

#===========================1. Stable temperatures (-5, 0, 5) + soil respiration only [non-linear mixed effects] END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#===========================2. Stable temperatures (-5, 0, 5) - root derived respiration
TrRoots <- TrLongerFull |> filter(Origin == "Roots" & StabTemp == 'Yes')
TrRoots$Locality <- as.factor(TrRoots$Locality)
# ---- Converting C0 from mg to ug as MonthInc from months to Time in days
TrRoots$C0 <- TrRoots$C0*1000 #ug per gram of dry soil
TrRoots$Time <- TrRoots$MonthsInc*30 #days
#checking the ranges
cbind(summary(TrRoots$RootResp), summary(TrRoots$C0))

# ---- Initial parameters guess for root decomposition ----
root_start_vals <- TrRoots %>%
  filter(RootResp > 0 & !is.na(RootResp)) |> 
  group_by(Temp) %>%
  group_modify(~{
    d <- .x
    fit <- lm(log(RootResp / C0) ~ Time, data = d)
    tibble(k = -coef(fit)["Time"])
  })

print(root_start_vals) ###Almost no temperature sensitivity - this must be tested explicitely first !!!!!!!!!!!!!!!!!!!!!!

# ---- Nested nonlinear mixed-effect models ----#
#1a. No temperature sensitivity at all
R_root_a <- gnls(
  RootResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = TrRoots,
  params = logk ~ 1,
  start = c(logk = log(0.0247)),
  control = gnlsControl(maxIter = 300, msMaxIter = 300),
  na.action = na.exclude
)
#1b. Exponential sensitivity
R_root_b <- update(R_root_a, params = logk ~ Temp,
                   start = c(logk = c(log(0.0247), 0)))
#1c. Gaussian relationship
R_root_c <- update(R_root_b, params = logk ~ Temp + I(Temp^2),
                   start = c(logk = c(log(0.0247), 0, 0)))

anova(R_root_a, R_root_b, R_root_c)
AIC(R_root_a, R_root_b, R_root_c)

summary(R_root_c) #The temperature sensitivity function has U shape !!!!!!!!!!!!!!!
Rt0 <- R_root_c
#2. Is the root decomposition different at different localities?
Rt1 <- update(Rt0, params = logk ~ Locality + Temp + I(Temp^2),
              start = c(logk = c(coef(Rt0)["logk.(Intercept)"], rep(0, 5),
                                 coef(Rt0)["logk.Temp"], coef(Rt0)["logk.I(Temp^2)"])))
#3. Is the temperature sensitivity of root respiration different at different localities?
Rt2 <- update(Rt1, params = logk ~ Locality + Locality:Temp + I(Temp^2),
              start = c(coef(Rt1)[1:6], rep(coef(Rt1)["logk.Temp"], 6), coef(Rt1)["logk.I(Temp^2)"]))
#4. Is the shape of temperature response different at different localities?
Rt3 <- update(Rt2, params = logk ~ Locality * Temp + Locality * I(Temp^2),
              start = c(coef(Rt2)[1:12], rep(coef(Rt2)[length(coef(Rt2))], 6)))

anova(Rt0, Rt1, Rt2, Rt3)
AIC(Rt0, Rt1, Rt2, Rt3)

# ---- fitting the heteroscedasticity on Rt0 ----
Rt0_diag <- tibble(fitted = fitted(Rt0), resid = residuals(Rt0),
                   Temp = TrRoots$Temp)

ggplot(Rt0_diag, aes(fitted, resid)) +
  geom_point(alpha=0.5) + geom_hline(yintercept=0, linetype="dashed") +
  labs(title = "Residuals vs fitted, Rt0")

ggplot(Rt0_diag, aes(factor(Temp), abs(resid))) + geom_boxplot() +
  labs(title = "Absolute residuals by Temp, Rt0")

TrRoots_cc <- TrRoots[complete.cases(TrRoots[, c("SoilResp", "Temp", "Time", "Locality")]), ]
nrow(TrRoots_cc)

Rt0_varIdent_temp <- gnls(
  RootResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = TrRoots_cc,
  params = logk ~ Temp + I(Temp^2),
  start = coef(Rt0),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)
anova(Rt0, Rt0_varIdent_temp)
AIC(Rt0, Rt0_varIdent_temp)
anova(Rt0_varIdent_temp)
R2M(Rt0_varIdent_temp)

#=======================alternative model specification
# Model B: single shared 'a', free but constant across all temperatures
# tests whether the a=1 assumption is supported AT ALL
RtB <- gnls(
  RootResp ~ exp(logk) * plogis(logit_a) * C0 * exp(-exp(logk) * Time),
  data = TrRoots_cc,
  params = list(logk ~ Temp + I(Temp^2), logit_a ~ 1),
  start = c(logk = coef(Rt0_varIdent_temp), logit_a = qlogis(0.8)),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)

# Model C: a differs specifically between -5C ("cold") and 0/5C ("warm")
# targeted test of your lysis-driven-availability hypothesis
TrRoots_cc$Cold <- factor(ifelse(TrRoots_cc$Temp == -5, "cold", "warm"))

RtC <- gnls(
  RootResp ~ exp(logk) * plogis(logit_a) * C0 * exp(-exp(logk) * Time),
  data = TrRoots_cc,
  params = list(logk ~ Temp + I(Temp^2), logit_a ~ Cold),
  start = c(logk = coef(RtB)[1:3], logit_a = c(coef(RtB)["logit_a"], 0)),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)

anova(Rt0_varIdent_temp, RtB)   # is a=1 rejected at all?
anova(RtB, RtC)                 # does a differ specifically at -5C?
logit_a_cold <- coef(RtC)["logit_a.(Intercept)"]
logit_a_warm <- coef(RtC)["logit_a.(Intercept)"] + coef(RtC)["logit_a.Coldwarm"]

a_cold <- plogis(logit_a_cold)
a_warm <- plogis(logit_a_warm)

c(a_cold = a_cold, a_warm = a_warm)
#============================================================================END!!!!!!!!!!!!!!!!!!!!!!!!!!

Rt0_boot_v2 <- boot_gnls(Rt0_varIdent_temp, nboot = 999)
boot_coefs_Rt0_v2 <- as.data.frame(Rt0_boot_v2$t)
names(boot_coefs_Rt0_v2) <- names(Rt0_boot_v2$t0)
boot_coefs_Rt0_v2 <- boot_coefs_Rt0_v2[complete.cases(boot_coefs_Rt0_v2), ]

b  <- boot_coefs_Rt0_v2[["logk.(Intercept)"]]
b1 <- boot_coefs_Rt0_v2[["logk.Temp"]]
b2 <- boot_coefs_Rt0_v2[["logk.I(Temp^2)"]]

T_min <- -b1 / (2*b2)          # vertex location
k_min <- exp(b - b1^2/(4*b2))  # rate at the vertex

tibble(
  metric = c("T_min", "k_min"),
  mean = c(mean(T_min), mean(k_min)),
  lwr  = c(quantile(T_min, 0.025), quantile(k_min, 0.025)),
  upr  = c(quantile(T_min, 0.975), quantile(k_min, 0.975))
)

k_roots_df <- lapply(temp_grid, function(T) {
  logk_roots  <- b1 + b1 * T + b2*T^2
  tibble(
    Temp = T,
    k_mean = c(mean(exp(logk_roots)), mean(exp(logk_roots))),
    k_lwr  = c(quantile(exp(logk_roots), 0.025), quantile(exp(logk_roots), 0.025)),
    k_upr  = c(quantile(exp(logk_roots), 0.975), quantile(exp(logk_roots), 0.975))
  )
}) %>% bind_rows()

print(k_roots_df)

ggplot(k_roots_df, aes(Temp, k_mean)) +
  geom_ribbon(aes(ymin = k_lwr, ymax = k_upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Temperature (°C)", y = "Kinetic constant k (1/d)") +
  mytheme + theme(legend.position = c(0.15, 0.85))

# ---- compare to T_cross ----
# T_min_draws should already exist from the earlier Rt0 extraction:
# T_min_draws <- -b1 / (2*b2)   using boot_coefs_Rt0_v2

T_diff <- T_cross - T_min

tibble(
  metric = c("T_cross (R2)", "T_min (Rt0)", "difference (cross - min)"),
  mean = c(mean(T_cross, na.rm=TRUE), mean(T_min, na.rm=TRUE), mean(T_diff, na.rm=TRUE)),
  lwr  = c(quantile(T_cross,0.025,na.rm=TRUE), quantile(T_min,0.025,na.rm=TRUE), quantile(T_diff,0.025,na.rm=TRUE)),
  upr  = c(quantile(T_cross,0.975,na.rm=TRUE), quantile(T_min,0.975,na.rm=TRUE), quantile(T_diff,0.975,na.rm=TRUE))
)

#===========================2. Stable temperatures (-5, 0, 5) - root derived respiration END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#===========================3. Freeze-thaw events
##=================3a. Control soils - applying M1 predictions to correct for possible difference in residual available organic carbon after 2 months incubation at -5 or -18
# M1_gnls was fit on Ctr — confirm complete cases
Ctr_cc <- Ctr[complete.cases(Ctr[, c("SoilResp", "Temp", "Time", "Locality")]), ]
nrow(Ctr_cc)

# ---- the predictions with bootstrap confidence interval are calculated first ----
M1_refit <- gnls(
  SoilResp ~ exp(logk) * C0 * exp(-exp(logk) * Time),
  data = Ctr_cc,
  params = list(C0 ~ Locality, logk ~ Temp),
  start = coef(M1_gnls),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)
# ---- Model coefficients, which will be used for the predictions ----
set.seed(123)
M1_boot <- boot_gnls(M1_refit, nboot = 999)
boot_coefs_M1 <- as.data.frame(M1_boot$t)
names(boot_coefs_M1) <- names(M1_boot$t0)
boot_coefs_M1_clean <- boot_coefs_M1[complete.cases(boot_coefs_M1), ]
nrow(boot_coefs_M1_clean)   # should be 999
# ---- Posthoc comparison of an active fraction of total organic carbon ----
localities <- levels(Ctr_cc$Locality)

C0_draws <- sapply(localities, function(loc) {
  if (loc == localities[1]) {
    boot_coefs_M1_clean[["C0.(Intercept)"]]
  } else {
    boot_coefs_M1_clean[["C0.(Intercept)"]] + boot_coefs_M1_clean[[paste0("C0.Locality", loc)]]
  }
})
colnames(C0_draws) <- localities

pairs_list <- combn(localities, 2, simplify = FALSE)

toc_by_soil <- Ctr_cc %>% distinct(Locality) %>%
  left_join(toc_mean_df, by = "Locality")   # from the earlier C0-vs-TOC plot exercise

toc_vec <- toc_by_soil$TOC_mean[match(localities, toc_by_soil$Locality)]

active_frac_draws <- sweep(C0_draws, 2, toc_vec, "/")   # C0/TOC per draw, per soil

frac_pairwise <- lapply(pairs_list, function(p) {
  diff_draws <- active_frac_draws[, p[1]] - active_frac_draws[, p[2]]
  tibble(
    comparison = paste(p[1], "-", p[2]),
    mean_diff = mean(diff_draws),
    lwr = quantile(diff_draws, 0.025),
    upr = quantile(diff_draws, 0.975),
    p_raw = 2 * min(mean(diff_draws > 0), mean(diff_draws < 0))
  )
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p_raw, method = "holm"))

print(frac_pairwise)

# ---- predictions calculated ----
localities <- levels(Ctr_cc$Locality)   
t2_days <- c(0, 30.44, 60.88)           # time since thaw: matches your 2/3/4 month measurements

C0_by_soil <- function(draws, locality) {
  base <- draws[["C0.(Intercept)"]]
  offset_name <- paste0("C0.Locality", locality)
  if (offset_name %in% names(draws)) {
    return(base + draws[[offset_name]])
  } else {
    return(base)   # this locality IS the reference level — no offset exists
  }
}

logk_at_temp <- function(draws, temp) {
  draws[["logk.(Intercept)"]] + draws[["logk.Temp"]] * temp
}

predict_FT_trajectory_v2 <- function(draws, locality, freeze_temp) {
  C0 <- C0_by_soil(draws, locality)
  k5 <- exp(logk_at_temp(draws, 5))
  
  if (freeze_temp == -5) {
    k_freeze <- exp(logk_at_temp(draws, -5))
  } else if (freeze_temp == -18) {
    k_freeze <- 0   # near-zero decay assumption
  }
  
  freeze_durations <- c(60.88, 30.44, 30.44)   # days: initial 2mo, then two 1mo refreezes
  thaw_duration <- 14
  
  C <- C0
  rates_at_measurement <- numeric(3)
  
  for (i in 1:3) {
    C <- C * exp(-k_freeze * freeze_durations[i])        # deplete during freeze phase
    rates_at_measurement[i] <- k5 * C * exp(-k5 * thaw_duration)  # rate at END of thaw window
    C <- C * exp(-k5 * thaw_duration)                     # deplete during the thaw window itself
  }
  
  rates_at_measurement
}

pred_all_v2 <- expand.grid(Locality = localities, Treatment = c("-5/5", "-18/5")) %>%
  rowwise() %>%
  group_split() %>%
  lapply(function(row) {
    freeze_temp <- ifelse(row$Treatment == "-5/5", -5, -18)
    sims <- t(apply(boot_coefs_M1_clean, 1, function(d)
      predict_FT_trajectory_v2(d, row$Locality, freeze_temp)))
    tibble(
      Locality = row$Locality,
      Treatment = row$Treatment,
      Time_since_thaw = c(0, 1, 2),   # measurement index, for plotting alignment with existing figure
      pred_mean = colMeans(sims, na.rm=TRUE),
      pred_lwr  = apply(sims, 2, quantile, 0.025, na.rm=TRUE),
      pred_upr  = apply(sims, 2, quantile, 0.975, na.rm=TRUE)
    )
  }) %>%
  bind_rows()

print(pred_all_v2)

# ---- Combine with measured data ----
FT_df <- TrLongerFull %>%
  filter(Origin == 'Control soil' & TempFactor %in% c("5", "-18/5", "-5/5")) %>%
  mutate(TempFactor = factor(TempFactor, levels = c("5", "-18/5", "-5/5")),
         Locality = factor(Locality))
FT_df$Time <- FT_df$MonthsInc*30

obs_FT <- FT_df %>%
  filter(TempFactor %in% c("-5/5", "-18/5")) %>%
  mutate(Time_since_thaw = round((Time - 60.88) / 30.44)) %>%   # months since thaw, from days
  rename(Treatment = TempFactor)

compare_df <- obs_FT %>%
  left_join(pred_all_v2, by = c("Locality", "Treatment", "Time_since_thaw")) %>%
  mutate(above_predicted = SoilResp > pred_upr,
         below_predicted = SoilResp < pred_lwr)

compare_df %>%
  group_by(Locality, Treatment) %>%
  summarize(n_above = sum(above_predicted, na.rm=TRUE),
            n_below = sum(below_predicted, na.rm=TRUE),
            n_total = n())

# ---- Plot the results ----
ggplot(pred_all, aes(Time_since_thaw, pred_mean)) +
  geom_ribbon(aes(ymin = pred_lwr, ymax = pred_upr), alpha = 0.2) +
  geom_line() +
  geom_point(data = obs_FT, aes(y = SoilResp), color = "red") +
  facet_grid(Locality ~ Treatment, scales = "free_y") +
  labs(x = "Months since thaw", y = "dR/dt",
       title = "Observed (red) vs. predicted (band) respiration after freeze-thaw")
# ---- Formal test of differences between predictions and observations ----
compare_df_v2 <- obs_FT %>%
  left_join(pred_all_v2, by = c("Locality", "Treatment", "Time_since_thaw")) %>%
  mutate(residual = SoilResp - pred_mean,
         above_band = SoilResp > pred_upr,
         below_band = SoilResp < pred_lwr)
# ---- The linear model fits the "residual" unexplained differences from expected trajectory defined in M1 model ----
resid_model_reduced <- lm(residual ~ Treatment + Time_since_thaw + Locality +
                            Treatment:Time_since_thaw + Treatment:Locality, data = compare_df_v2)
anova(resid_model_reduced)
# ---- Visualize ----
compare_df_v2 %>%
  group_by(Treatment, Time_since_thaw) %>%
  summarize(mean_resid = mean(residual, na.rm=TRUE),
            se_resid = sd(residual, na.rm=TRUE)/sqrt(n()),
            n = n(),
            .groups = "drop")

compare_df_v2 %>%
  group_by(Treatment, Time_since_thaw) %>%
  summarize(mean_resid = mean(residual, na.rm=TRUE),
            se_resid = sd(residual, na.rm=TRUE)/sqrt(n()),
            .groups = "drop") %>%
  ggplot(aes(Time_since_thaw + 2, mean_resid, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", lwd = 1) +
  geom_pointrange(aes(ymin = mean_resid - se_resid, ymax = mean_resid + se_resid),
                  position = position_dodge(width = 0.1), size = 0.8) +
  geom_line(position = position_dodge(width = 0.1)) +
  labs(y = "Observed − Predicted (dR/dt)", x = "Time of incubation (months)") +
  theme_min + theme(legend.position = c(0.8, 0.8))
##=================3a. Control soils - applying M1 predictions END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##=================3b. Soil respiration affected by root addition
# ---- R2 model predictions ----
C0_by_soil_R2 <- function(draws, locality) {
  base <- draws[["C0.(Intercept)"]]
  offset_name <- paste0("C0.Locality", locality)
  if (offset_name %in% names(draws)) base + draws[[offset_name]] else base
}

logk_at_temp_R2_rootpresent <- function(draws, temp) {
  draws[["logk.(Intercept)"]] +
    draws[["logk.Temp"]] * temp +
    draws[["logk.OriginSoil with roots"]] +
    draws[["logk.Temp:OriginSoil with roots"]] * temp
}

predict_FT_trajectory_R2 <- function(draws, locality, freeze_temp) {
  C0 <- C0_by_soil_R2(draws, locality)
  k5 <- exp(logk_at_temp_R2_rootpresent(draws, 5))
  k_freeze <- if (freeze_temp == -5) exp(logk_at_temp_R2_rootpresent(draws, -5)) else 0
  
  freeze_durations <- c(60.88, 30.44, 30.44)
  thaw_duration <- 14
  
  C <- C0
  rates <- numeric(3)
  for (i in 1:3) {
    C <- C * exp(-k_freeze * freeze_durations[i])
    rates[i] <- k5 * C * exp(-k5 * thaw_duration)
    C <- C * exp(-k5 * thaw_duration)
  }
  rates
}

# ---- Comparison ----
RootFT_soil_pred <- expand.grid(Locality = localities, Treatment = c("-5/5","-18/5")) %>%
  rowwise() %>% group_split() %>%
  lapply(function(row) {
    ft <- ifelse(row$Treatment == "-5/5", -5, -18)
    sims <- t(apply(boot_coefs_v2, 1, function(d) predict_FT_trajectory_R2(d, row$Locality, ft)))
    tibble(Locality = row$Locality, Treatment = row$Treatment, Time_since_thaw = c(0,1,2),
           pred_mean = colMeans(sims, na.rm=TRUE),
           pred_lwr = apply(sims, 2, quantile, 0.025, na.rm=TRUE),
           pred_upr = apply(sims, 2, quantile, 0.975, na.rm=TRUE))
  }) %>% bind_rows()

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

compare_df_R2_FT <- obs_FT_rootadded %>%
  left_join(RootFT_soil_pred, by = c("Locality","Treatment","Time_since_thaw")) %>%
  mutate(residual = SoilResp - pred_mean)

resid_model_R2 <- lm(residual ~ Treatment * Time_since_thaw + Locality, data = compare_df_R2_FT)
anova(resid_model_R2)

# Predictions are biased due to overfit of R2, alternative approach is selected:

# ---- Alternative direct comparison of respiration rates ----
obs_compare <- bind_rows(
  obs_FT %>% mutate(Root = "absent"),
  obs_FT_rootadded %>% mutate(Root = "present")
) %>%
  mutate(Root = factor(Root, levels = c("absent","present")),
         Treatment = factor(Treatment, levels = c("-5/5","-18/5")),
         Locality = factor(Locality))

# ---- statistical models ----
resid_model_full <- lm(SoilResp ~ Root * Treatment * Time_since_thaw, data = obs_compare)
anova(resid_model_full)
resid_model_reduced <- lm(SoilResp ~ Root + Treatment + Time_since_thaw +
                            Root:Treatment + Root:Time_since_thaw + Treatment:Time_since_thaw,
                          data = obs_compare)
anova(resid_model_reduced)
anova(resid_model_reduced, resid_model_full)   # is the 3-way interaction worth keeping?

# ---- visualize ----
obs_compare %>%
  group_by(Root, Treatment, Time_since_thaw) %>%
  summarize(mean_resp = mean(SoilResp, na.rm=TRUE),
            se_resp = sd(SoilResp, na.rm=TRUE)/sqrt(n()), .groups="drop")

obs_compare %>%
  group_by(Root, Treatment, Time_since_thaw) %>%
  summarize(mean_resp = mean(SoilResp, na.rm=TRUE),
            se_resp = sd(SoilResp, na.rm=TRUE)/sqrt(n()), .groups="drop") %>%
  ggplot(aes(Time_since_thaw + 2, mean_resp, color = Root, linetype = Treatment)) +
  geom_pointrange(aes(ymin = mean_resp - se_resp, ymax = mean_resp + se_resp),
                  position = position_dodge(width=0.1)) +
  geom_line(position = position_dodge(width=0.1)) +
  labs(y = "Respiration", x = "Time of incubation (months)") +
  mytheme + theme(legend.direction = 'horizontal', legend.position = c(0.7, 0.7))


##=================3b. Soil respiration affected by root addition END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##=================3c. Root derived respiration
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

ggplot(RootFT_obs %>% group_by(Treatment, Time_since_thaw) %>%
         summarize(mean_resp = mean(RootResp, na.rm=TRUE),
                   se_resp = sd(RootResp, na.rm=TRUE)/sqrt(n()), .groups="drop"),
       aes(Time_since_thaw+2, mean_resp, color = Treatment)) +
  geom_pointrange(aes(ymin = mean_resp-se_resp, ymax = mean_resp+se_resp)) +
  geom_line() +
  labs(y = "Root derived respiration", x = "Time of incubation (months)") +
  mytheme + theme(legend.direction = 'horizontal', legend.position = c(0.7, 0.7))

resid_root_model <- lm(RootResp ~ Treatment * Time_since_thaw, data = RootFT_obs)
anova(resid_root_model)

month0_test <- obs_compare %>% filter(Time_since_thaw == 0, Treatment == "-18/5")
t.test(SoilResp ~ Root, data = month0_test)

# ---- Testing the -5 anomaly against freeze-thaw treatments - the anomaly is caused by heavy frost lysis ----
# ---- comparing -5 with -5/5 ----
k_neg5 <- exp(coef(Rt0_varIdent_temp)["logk.(Intercept)"] + coef(Rt0_varIdent_temp)["logk.Temp"]*(-5) +
                coef(Rt0_varIdent_temp)["logk.I(Temp^2)"]*25)
k_5    <- exp(coef(Rt0_varIdent_temp)["logk.(Intercept)"] + coef(Rt0_varIdent_temp)["logk.Temp"]*5 +
                coef(Rt0_varIdent_temp)["logk.I(Temp^2)"]*25)

kinetic_ratio <- (k_5 / k_neg5) * exp(-k_5 * 14)
obs_stable_neg5_m2 <- obs_stable_neg5$obs_mean[obs_stable_neg5$Time == 60]

kinetic_only_pred_5_5 <- kinetic_ratio * obs_stable_neg5_m2

obs_5_5_m0 <- RootFT_obs %>% filter(Treatment == "-5/5", Time_since_thaw == 0) %>%
  summarize(mean_resp = mean(RootResp, na.rm=TRUE)) %>% pull()

cat("Kinetic-only prediction for -5/5 at month 2:", kinetic_only_pred_5_5, "\n")
cat("Observed -5/5 at month 2:", obs_5_5_m0, "\n")
cat("Excess ratio (observed / kinetic-only):", obs_5_5_m0 / kinetic_only_pred_5_5, "\n")
