####################################################################################################
###                        Tatra Mountains soil cores incubation experiment                      ###
######                          Statistical analysis alternative                             #######
####################################################################################################
#==================Statistical tests performed previously was identified as biased after carefull
#discussion. Since it is not possible to patch the previous analyses, the fresh start is conducted
#data frames and statistical models are re-used here to avoid extra work.

####################################################################################################
#========================================Statistical analysis======================================#
####################################################################################################
#===========================1. Stable temperatures (-5, 0, 5) + soil respiration only [non-linear mixed effects]
##=================1a. Control soils - nested models
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

# ---- Does the temperature affect a? It also should because the freeze-thaw inflates respiration ----
N2 <- update(N1, params = list(logk ~ Temp, logit_a ~ Locality + Temp),
             start = c(coef(N1), 0))

summary(N2)
anova(N1, N2)
AIC(N1, N2)
R2M(N2)
anova(N2, type = 'sequential', adjustSigma = T)

# ---- The temperature sensitivity looks off being negative, let's confirm that the there is not systematic bias visually ----
Ctr_cc$fitted_N2 <- fitted(N2)

ggplot(Ctr_cc, aes(Time, SoilResp, color=factor(Temp))) +
  geom_point(alpha=0.5) +
  geom_line(aes(y=fitted_N2, group=interaction(Locality,Temp))) +
  facet_wrap(~Locality, scales="free_y")

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
##=================Figure - measured vs predicted respiration rates
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

RootFreeStable %>%
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

#count number of measurements within the confidence interval of predictions
Ctr_cc$In <- numeric(length = nrow(Ctr_cc))
for(i in 1:nrow(Ctr_cc)){
  #Model parameters for row i
  loc_i <- Ctr_cc$Locality[i]
  k_i = exp(boot_coefs_N2[['logk.(Intercept)']] + boot_coefs_N2[['logk.Temp']]*Ctr_cc$Temp[i])
  if(loc_i == 'FU'){
    alfa_i =  plogis(boot_coefs_N2[['logit_a.(Intercept)']] + 
                       boot_coefs_N2[['logit_a.Temp']]*Ctr_cc$Temp[i])
  } else {
    alfa_i = plogis(boot_coefs_N2[['logit_a.(Intercept)']] + 
                      boot_coefs_N2[[paste0('logit_a.Locality', loc_i)]] + 
                      boot_coefs_N2[['logit_a.Temp']]*Ctr_cc$Temp[i])
  }
  
  #Prediction for row i
  interval_i <- as.numeric(quantile(x = Ctr_cc$TOC_anchor[i]*
                                      k_i*alfa_i*exp(-k_i*Ctr_cc$Time[i]),
                                    probs = c(0.025, 0.975)))
  
  Ctr_cc$In[i] <- ifelse(dplyr::between(x = Ctr_cc$SoilResp[i], left = interval_i[1], right = interval_i[2]),
                         1, 0)
  
}  
sum(Ctr_cc$In)/nrow(Ctr_cc)
#=============Residual fraction of initial activity
for(i in c(-5, 0, 5)){
  print(median(exp(-exp(boot_coefs_N2[['logk.(Intercept)']] + boot_coefs_N2[['logk.Temp']]*i)*120)*100))
}

##=================1a. Control soils - nested models END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##=================1b. Control soils and soils with roots - refiting original models
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

# ---- The best models is surprisingly Rc3 - summarize and bootstrap ----
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

TruePriming <- data.frame(Legend = rep(c("Control", "Root-added"), each = nrow(boot_coefs_Rc3)),
                          k = c(exp(boot_coefs_Rc3[, 'logk.(Intercept)']),
                                exp(boot_coefs_Rc3[, 'logk.(Intercept)'] + 
                                      boot_coefs_Rc3[, 'logk.OriginSoil with roots']))) 

library(gridExtra)

grid.arrange(
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

##=================Figure - measured vs predicted respiration rates
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

RootAmendedStable %>%
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

#count number of measurements within the confidence interval of predictions
TrSoil_cc2$In <- numeric(length = nrow(TrSoil_cc2))
for(i in 1:nrow(TrSoil_cc2)){
  #Model parameters for row i
  loc_i <- TrSoil_cc2$Locality[i]
  k_i = exp(boot_coefs_Rc3[['logk.(Intercept)']] + boot_coefs_Rc3[['logk.OriginSoil with roots']] + 
              (boot_coefs_Rc3[['logk.Temp']] +
                 boot_coefs_Rc3[['logk.Temp:OriginSoil with roots']])*TrSoil_cc2$Temp[i])
  if(loc_i == 'FU'){
    alfa_i =  plogis(boot_coefs_Rc3[['logit_a.(Intercept)']] + boot_coefs_Rc3[['logit_a.OriginSoil with roots']] +
                       boot_coefs_Rc3[['logit_a.Temp']]*TrSoil_cc2$Temp[i])
  } else {
    alfa_i = plogis(boot_coefs_Rc3[['logit_a.(Intercept)']] + boot_coefs_Rc3[['logit_a.OriginSoil with roots']] +
                      boot_coefs_Rc3[[paste0('logit_a.Locality', loc_i)]] + 
                      boot_coefs_Rc3[['logit_a.Temp']]*TrSoil_cc2$Temp[i])
  }
  
  
  #Prediction for row i
  interval_i <- as.numeric(quantile(x = TrSoil_cc2$TOC_anchor[i]*
                                      k_i*alfa_i*exp(-k_i*TrSoil_cc2$Time[i]),
                                    probs = c(0.025, 0.975)))
  
  TrSoil_cc2$In[i] <- ifelse(dplyr::between(x = TrSoil_cc2$SoilResp[i], left = interval_i[1], right = interval_i[2]),
                         1, 0)
  
}  
sum(TrSoil_cc2[TrSoil_cc2$Origin == "Soil with roots", "In"])/nrow(TrSoil_cc2[TrSoil_cc2$Origin == "Soil with roots", ])
#=============Residual fraction of initial activity
for(i in c(-5, 0, 5)){
  print(median(exp(-exp(boot_coefs_Rc3[['logk.(Intercept)']] + boot_coefs_Rc3[['logk.OriginSoil with roots']] + 
                          (boot_coefs_Rc3[['logk.Temp']] +
                             boot_coefs_Rc3[['logk.Temp:OriginSoil with roots']])*i)*120)*100))
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

#===========================1. Stable temperatures (-5, 0, 5) + soil respiration only [non-linear mixed effects] END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#===========================2. Stable temperatures (-5, 0, 5) - root derived respiration
# ---- Figure showing fraction of root-derived respiration for supplementary material ----
TrRoots_cc %>%
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

# ---- the new null model with available fraction of root material directly acknowledged ---- 
Pa_k0 <- gnls(
  RootResp ~ exp(logk) * plogis(logit_a) * C0 * exp(-exp(logk) * Time),
  data = TrRoots_cc,
  params = list(logk ~ 1, logit_a ~ 1),
  start = c(logk = log(0.05), logit_a = qlogis(0.15)),
  weights = varIdent(form = ~ 1 | Temp),
  control = gnlsControl(maxIter = 300, msMaxIter = 300)
)
# ---- Expoentil temperature function ----
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

# ---- the best model is  Pc_k_loc - summarize and bootstrap ----
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

##=================Figure - measured vs predicted respiration rates
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

RootStable %>%
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

#count number of measurements within the confidence interval of predictions
TrRoots_cc$In <- numeric(length = nrow(TrRoots_cc))
for(i in 1:nrow(TrRoots_cc)){
  #Model parameters for row i
  loc_i <- TrRoots_cc$Locality[i]
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
  interval_i <- as.numeric(quantile(x = TrRoots_cc$C0[i]*
                                      k_i*alfa_i*exp(-k_i*TrRoots_cc$Time[i]),
                                    probs = c(0.025, 0.975), na.rm = T))
  
  TrRoots_cc$In[i] <- ifelse(dplyr::between(x = TrRoots_cc$SoilResp[i], left = interval_i[1], right = interval_i[2]),
                         1, 0)
  
}  
sum(TrRoots_cc$In, na.rm = T)/nrow(TrRoots_cc[!is.na(TrRoots_cc$In), ])*100
#=============Residual fraction of initial activity
for(i in c(-5, 0, 5)){
  print(mean(exp(-exp(boot_coefs_Pc_k[['logk.(Intercept)']] + boot_coefs_Pc_k[['logk.Temp']]*i)*60)*100))
}

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


#===========================2. Stable temperatures (-5, 0, 5) - root derived respiration END!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#===========================3. Freeze-thaw events
##=================3a. Control soils - applying N2 predictions to correct for possible difference in residual available organic carbon after 2 months incubation at -5 or -18
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

# ---- let's do the comparison against stable -5 and 5C as well ----

anova(aov(SoilResp ~ TempFactor*Locality, data = TrLongerFull, 
          subset = (MonthsInc == 2 & Origin == 'Soil' & TempFactor != '-18/5' & Temp != '0')))

# ---- Time 0 (2 months) fits the prediction cleanly. this is our well defined reference point ----
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


# ---- supplementary figure ----
trajDFrame$Series <- c('-5/5')
trajDFrame18$Series <- c('-18/5')

rbind(trajDFrame, trajDFrame18) |> group_by(Series, Scenario, Locality, Month) |> 
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
  theme(axis.title.x = element_blank()) +
  labs(y = c("Relative effect of alternating thermal regime (nondimensional)"),
       color = c("Time (months)"), fill = c("Time (months)"))

# ---- At the end, we could compare -5/5 and -18/5 treatments directly ----
ft_raw <- obs_FT %>% filter(Treatment %in% c("-5/5","-18/5"))

# overall test across all three months
ft_raw_model <- lm(SoilResp ~ Treatment * Time_since_thaw*Locality, data = ft_raw, subset = Locality != "RA")
anova(ft_raw_model)
summary(ft_raw_model)

# direct, targeted comparison at month 2 specifically
t.test(SoilResp ~ Treatment, data = ft_raw %>% filter(Time_since_thaw == 0))

# and the full by-month summary for context
ft_raw %>% group_by(Treatment, Time_since_thaw) %>%
  summarize(mean_resp = mean(SoilResp, na.rm=TRUE), se = sd(SoilResp,na.rm=TRUE)/sqrt(n()), .groups="drop")

# ---- Figure ----
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

library(ggnewscale)

RootFreeRef <- Ctr_cc |> group_by(Locality, Temp, MonthsInc) |> 
  summarise(yref = mean(SoilResp, na.rm = T), yrefSD = sd(SoilResp, na.rm = T)) |> 
  filter(Temp != 0) |> filter(MonthsInc != 3) |>  
  mutate(TempF = factor(Temp, levels = c(-5, 5), labels = c("-5˚C", '+5˚C')))

pred_month2_FT2 <- rbind(pred_month2_FT, pred_month2_FT)
pred_month2_FT2$MonthsInc <- rep(c(2, 4), each = nrow(pred_month2_FT)) 
pred_month2_FT2[(nrow(pred_month2_FT)+1):nrow(pred_month2_FT2), c(2:4)] <- NA

TrLongerFull |> filter(Roots == 'No' & StabTemp == 'No' & MonthsInc != 3) |> 
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
  theme(axis.title.x = element_blank(), legend.position = c(0.1, 0.9)) +
  scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                              expression(atop("2160 m", "meadow")),
                              expression(atop("2130 m", "meadow")),
                              expression(atop("1780 m", "meadow")),
                              expression(atop("1470 m", "forest")),
                              expression(atop("1330 m", "forest"))))


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

##=================3a. Control soils - applying N2 predictions END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##=================3b. Soil respiration affected by root addition
# ---- The same pipeline is implemented ----
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

# ---- Defining and testing the excess ratio using wilcoxen test ----
compare_5_5_Rc3 <- pred_5_5_val %>%
  left_join(obs_5_5_Rc3, by = c("Locality","Root")) %>%
  mutate(within_band = obs_mean >= pred_lwr & obs_mean <= pred_upr,
         excess_ratio = obs_mean / pred_mean)

print(compare_5_5_Rc3, n=30)

# ---- the formal test  - three different to be absolutely sure ----
excess_wide <- compare_5_5_Rc3 %>%
  select(Locality, Root, excess_ratio) %>%
  tidyr::pivot_wider(names_from = Root, values_from = excess_ratio)

print(excess_wide)

t.test(excess_wide$present, excess_wide$absent, paired = TRUE)
wilcox.test(excess_wide$present, excess_wide$absent, paired = TRUE)

# ---- re-confirming with SC 0, which makes the ratio undefined ----
excess_wide_log <- excess_wide %>%
  mutate(log_present = log(present + 0.01), log_absent = log(absent + 0.01))

t.test(excess_wide_log$log_present, excess_wide_log$log_absent, paired = TRUE)

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

compare_18_5_Rc3 <- pred_18_5_val %>%
  left_join(obs_18_5_Rc3, by=c("Locality","Root")) %>%
  mutate(excess_ratio = obs_mean / pred_mean)

print(compare_18_5_Rc3, n=30)

# ---- comparison ----
excess_compare <- excess_wide %>%
  rename(absent_5_5 = absent, present_5_5 = present) %>%
  left_join(
    compare_18_5_Rc3 %>% select(Locality, Root, excess_ratio) %>%
      tidyr::pivot_wider(names_from=Root, values_from=excess_ratio) %>%
      rename(absent_18_5 = absent, present_18_5 = present),
    by="Locality"
  )

print(excess_compare)

t.test(excess_compare$absent_18_5, excess_compare$absent_5_5, paired=TRUE)
t.test(excess_compare$present_18_5, excess_compare$present_5_5, paired=TRUE)

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

# # combine with the control series already built (rel_decline_full has -5/5, -18/5, stable_-5, stable_5)
# rel_decline_all <- bind_rows(
#   rel_decline_full %>% filter(Series %in% c("-5/5","-18/5")) %>% mutate(Series = paste0(Series, "_control")),
#   rel_decline_root
# )
# 
# rel_decline_all %>% count(Series, Time_since_thaw)
# 
# rel_decline_all %>% group_by(Series, Time_since_thaw) %>%
#   summarize(mean_rel = mean(rel_change, na.rm=TRUE), se = sd(rel_change,na.rm=TRUE)/sqrt(n()), .groups="drop") %>%
#   tidyr::pivot_wider(names_from=Time_since_thaw, values_from=c(mean_rel,se))

# ---- formal test ----
# rel_decline_all <- rel_decline_all %>%
#   mutate(Root = ifelse(grepl("root", Series), "present", "absent"),
#          Treatment = ifelse(grepl("18", Series), "-18/5", "-5/5"))

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


# rel_decline_root_clean %>% filter(Series=="-5/5_root") %>%
#   group_by(Time_since_thaw) %>%
#   summarize(mean_rel = mean(rel_change, na.rm=TRUE), se = sd(rel_change,na.rm=TRUE)/sqrt(n()), n=n(), .groups="drop")
# 
# rel_decline_all_clean <- bind_rows(
#   rel_decline_full %>% filter(Series %in% c("-5/5","-18/5")) %>% mutate(Series = paste0(Series, "_control")),
#   rel_decline_root_clean
# ) %>%
#   mutate(Root = ifelse(grepl("root", Series), "present", "absent"),
#          Treatment = ifelse(grepl("18", Series), "-18/5", "-5/5"))
# 
# rel_decline_all_clean %>% group_by(Series, Time_since_thaw) %>%
#   summarize(mean_rel = mean(rel_change, na.rm=TRUE), se=sd(rel_change,na.rm=TRUE)/sqrt(n()), .groups="drop") %>%
#   tidyr::pivot_wider(names_from=Time_since_thaw, values_from=c(mean_rel,se))
# 
# rel_time_model_clean <- lm(rel_change ~ Treatment * Time_since_thaw,
#                            data = rel_decline_all_clean %>% filter(Time_since_thaw > 0 & Roots == "Yes"))
# 
# anova(rel_time_model_clean)
# summary(rel_time_model_clean)

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

# ---- supplementary figure ----
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

rbind(trajDFrame_soc, trajDFrame_soc18) |> group_by(Series, Scenario, Locality, Month, .drop = F) |> 
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
  theme(axis.title.x = element_blank()) +
  labs(y = c("Relative effect of alternating thermal regime (nondimensional)"),
       color = c("Time (months)"), fill = c("Time (months)"))




# ---- At the end, we could compare -5/5 and -18/5 treatments directly ----
ft_raw <- obs_FT %>% filter(Treatment %in% c("-5/5","-18/5"))

# overall test across all three months
ft_raw_model <- lm(SoilResp ~ Treatment * Time_since_thaw*Locality, data = ft_raw, subset = Locality != "RA")
anova(ft_raw_model)
summary(ft_raw_model)

# direct, targeted comparison at month 2 specifically
t.test(SoilResp ~ Treatment, data = ft_raw %>% filter(Time_since_thaw == 0))

# and the full by-month summary for context
ft_raw %>% group_by(Treatment, Time_since_thaw) %>%
  summarize(mean_resp = mean(SoilResp, na.rm = T))


# ---- figure: -5/5 and -18/5 in different facets for clarity ----
obs_FT_all <- TrLongerFull |> filter(StabTemp == 'No' & MonthsInc == 2) |> 
  group_by(Locality, Roots, TempFactor) |> summarise(y = median(SoilResp, na.rm = T),
                                                     ySD = sd(SoilResp, na.rm = T)) |> 
  mutate(Root = recode(Roots, Yes = 'present', No = 'absent'))

anova(lm(SoilResp ~ Roots, TrLongerFull, subset = StabTemp == 'No' & MonthsInc == 2))
summary(lm(SoilResp ~ Roots, TrLongerFull, subset = StabTemp == 'No' & MonthsInc == 2))

pred_FT_all <- bind_rows(pred_5_5_val, pred_18_5_val) |> mutate(TempFactor = rep(c('-5/5', '-18/5'), each = 12))
print(pred_FT_all, n = 24)

bind_rows(obs_FT_all, pred_FT_all) |> 
  ggplot(aes(factor(Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA")), y)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 1), aes(fill = Root),
           color = 'black') +
  mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD, color = Root), 
                          position = position_dodge(width = 1), show.legend = F) +
  scale_fill_manual(values = c('grey30', 'firebrick'), labels = c("Control", "Root-added")) +
  scale_color_manual(values = c('black', 'black')) +
  labs(color = 'Roots', fill = "Roots",
       y = expression(paste("Respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")"))) +
  new_scale_fill() +
  new_scale_color() +
  facet_wrap(.~factor(TempFactor, levels = c('-5/5', '-18/5')), scales = 'free_y') +
  geom_point(cex = 4, pch = 21, aes(y = pred_mean, fill = Root), position = position_dodge(width = 1),
             show.legend = F) + 
  geom_errorbar(aes(ymin = pred_lwr, ymax = pred_upr, colour = Root), 
                position = position_dodge(width = 1), width = 0.1, show.legend = F) +
  scale_fill_manual(values = c('red', 'red')) +
  scale_color_manual(values = c('red', 'red')) +
  theme(axis.title.x = element_blank(), legend.position = c(0.1, 0.85), legend.title = element_blank()) +
  scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                              expression(atop("2160 m", "meadow")),
                              expression(atop("2130 m", "meadow")),
                              expression(atop("1780 m", "meadow")),
                              expression(atop("1470 m", "forest")),
                              expression(atop("1330 m", "forest"))))


RootAmendedRef <- TrSoil_cc |>  filter(Origin == "Soil with roots") |> group_by(Locality, Temp, MonthsInc) |> 
  summarise(yref = mean(SoilResp, na.rm = T), yrefSD = sd(SoilResp, na.rm = T)) |> 
  filter(Temp != 0) |> filter(MonthsInc != 3) |> 
  mutate(TempF = factor(Temp, levels = c(-5, 5), labels = c("-5˚C", '+5˚C')))

pred_FT_all_2 <- rbind(pred_FT_all[pred_FT_all$Root == 'present', ], pred_FT_all[pred_FT_all$Root == 'present', ])
pred_FT_all_2 <- pred_FT_all_2[, -2]
pred_FT_all_2$MonthsInc <- rep(c(2, 4), each = 12) 
pred_FT_all_2[13:nrow(pred_FT_all_2), c(2:4)] <- NA

TrLongerFull |> filter(Origin == 'Soil with roots' & StabTemp == 'No' & MonthsInc != 3) |> 
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
  theme(axis.title.x = element_blank(), legend.position = c(0.1, 0.9)) +
  scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                              expression(atop("2160 m", "meadow")),
                              expression(atop("2130 m", "meadow")),
                              expression(atop("1780 m", "meadow")),
                              expression(atop("1470 m", "forest")),
                              expression(atop("1330 m", "forest"))))

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

##=================3b. Soil respiration affected by root addition END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##=================3c. Root derived respiration
# ---- Figure showing fraction of root-derived respiration for supplementary material ----
TrLongerFull %>% filter(Origin == 'Roots' & StabTemp == 'No') |> 
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
  theme(legend.position = c(0.08, 0.36))
# ---- The same pipeline again ----
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

# ---- formal test ----
compare_5_5_root <- pred_5_5_root_val %>%
  left_join(obs_FT_root %>% filter(Treatment=="-5/5"), by="Locality") %>%
  mutate(excess_ratio = obs_mean / pred_mean)

compare_18_5_root <- pred_18_5_root_val %>%
  left_join(obs_FT_root %>% filter(Treatment=="-18/5"), by="Locality") %>%
  mutate(excess_ratio = obs_mean / pred_mean)

print(compare_5_5_root)
print(compare_18_5_root)

t.test(log(compare_5_5_root$excess_ratio), mu = 0)   # log scale, since ratios are naturally right-skewed
wilcox.test(compare_5_5_root$excess_ratio, mu = 1)

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

# ---- Figure ----
pred_month2_FT_root <- bind_rows(pred_5_5_root_val, pred_18_5_root_val) |> 
  mutate(TempFactor = rep(c("-5/5", "-18/5"), each = 6)) 

TrLongerFull |> filter(Origin == 'Roots' & StabTemp == 'No' & MonthsInc == 2) |> 
  group_by(Locality, TempFactor) |> summarise(y = median(RootResp, na.rm = T),
                                              ySD = sd(RootResp, na.rm = T)) |> 
  inner_join(pred_month2_FT_root, by = c('Locality', "TempFactor")) |> 
  ggplot(aes(factor(Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA")), y)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 1), aes(fill = TempFactor),
           color = 'black') +
  mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD, color = TempFactor), 
                          position = position_dodge(width = 1), show.legend = F) +
  scale_fill_manual(values = c('white', 'grey')) +
  scale_color_manual(values = c('black', 'black')) +
  labs(color = 'Temperature regime', fill = "Temperature regime",
       y = expression(paste("Respiration (", mu, "g C-C", O[2], " g ", DW^{-1}, "da", y^{-1}, ")")),
       title = "B)") +
  new_scale_fill() +
  new_scale_color() +
  geom_point(cex = 4, pch = 21, aes(y = pred_mean, fill = TempFactor), position = position_dodge(width = 1),
             show.legend = F) + 
  geom_errorbar(aes(ymin = pred_lwr, ymax = pred_upr, colour = TempFactor), 
                position = position_dodge(width = 1), width = 0.1, show.legend = F) +
  scale_fill_manual(values = c('red', 'red')) +
  scale_color_manual(values = c('red', 'red')) +
  theme(axis.title.x = element_blank(), legend.position = c(0.15, 0.8)) +
  scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                              expression(atop("2160 m", "meadow")),
                              expression(atop("2130 m", "meadow")),
                              expression(atop("1780 m", "meadow")),
                              expression(atop("1470 m", "forest")),
                              expression(atop("1330 m", "forest"))))

grid.arrange(TrLongerFull |> filter(Origin == 'Roots' & StabTemp == 'No' & MonthsInc == 2) |> 
               group_by(Locality, TempFactor) |> summarise(y = median(SoilResp, na.rm = T),
                                                           ySD = sd(SoilResp, na.rm = T)) |> 
               inner_join(pred_month2_FT, by = c('Locality', "TempFactor")) |> 
               ggplot(aes(factor(Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA")), y)) + 
               geom_bar(stat = 'identity', position = position_dodge(width = 1), aes(fill = TempFactor),
                        color = 'black', show.legend = F) +
               mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD, color = TempFactor), 
                                       position = position_dodge(width = 1), show.legend = F) +
               scale_fill_manual(values = c('white', 'grey')) +
               scale_color_manual(values = c('black', 'black')) +
               labs(color = 'Temperature regime', fill = "Temperature regime",
                    y = expression(paste("Respiration (", mu, "g C-C", O[2], " g ", DW^{-1}~d^{-1}, ")")),
                    title = "A)") +
               new_scale_fill() +
               new_scale_color() +
               geom_point(cex = 4, pch = 21, aes(y = pred_mean, fill = TempFactor), position = position_dodge(width = 1),
                          show.legend = F) + 
               geom_errorbar(aes(ymin = pred_lwr, ymax = pred_upr, colour = TempFactor), 
                             position = position_dodge(width = 1), width = 0.1, show.legend = F) +
               scale_fill_manual(values = c('red', 'red')) +
               scale_color_manual(values = c('red', 'red')) +
               theme(axis.title.x = element_blank(), legend.position = c(0.15, 0.8), axis.text.x = element_blank()) +
               scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                                           expression(atop("2160 m", "meadow")),
                                           expression(atop("2130 m", "meadow")),
                                           expression(atop("1780 m", "meadow")),
                                           expression(atop("1470 m", "forest")),
                                           expression(atop("1330 m", "forest")))),
             TrLongerFull |> filter(Origin == 'Roots' & StabTemp == 'No' & MonthsInc == 2) |> 
               group_by(Locality, TempFactor) |> summarise(y = median(RootResp, na.rm = T),
                                                           ySD = sd(RootResp, na.rm = T)) |> 
               inner_join(pred_month2_FT_root, by = c('Locality', "TempFactor")) |> 
               ggplot(aes(factor(Locality, levels = c("SC", "VW", "NW", "FU", "JA", "RA")), y)) + 
               geom_bar(stat = 'identity', position = position_dodge(width = 1), aes(fill = TempFactor),
                        color = 'black') +
               mytheme + geom_errorbar(aes(ymin = y, ymax = y + ySD, color = TempFactor), 
                                       position = position_dodge(width = 1), show.legend = F) +
               scale_fill_manual(values = c('white', 'grey')) +
               scale_color_manual(values = c('black', 'black')) +
               labs(color = 'Temperature regime', fill = "Temperature regime",
                    y = expression(paste("Respiration (", mu, "g C-C", O[2], " g ", DW^{-1}~d^{-1}, ")")),
                    title = "B)") +
               new_scale_fill() +
               new_scale_color() +
               geom_point(cex = 4, pch = 21, aes(y = pred_mean, fill = TempFactor), position = position_dodge(width = 1),
                          show.legend = F) + 
               geom_errorbar(aes(ymin = pred_lwr, ymax = pred_upr, colour = TempFactor), 
                             position = position_dodge(width = 1), width = 0.1, show.legend = F) +
               scale_fill_manual(values = c('red', 'red')) +
               scale_color_manual(values = c('red', 'red')) +
               theme(axis.title.x = element_blank(), legend.position = c(0.15, 0.8)) +
               scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                                           expression(atop("2160 m", "meadow")),
                                           expression(atop("2130 m", "meadow")),
                                           expression(atop("1780 m", "meadow")),
                                           expression(atop("1470 m", "forest")),
                                           expression(atop("1330 m", "forest")))),
             ncol = 1)

# ---- temporal trajectory check ----
library(dplyr)

rel_decline_root_derived <- RootFT_obs %>%
  filter(TempFactor == "-5/5" | TempFactor == "-18/5") %>%
  group_by(Locality, TempFactor) %>%
  mutate(month2_val = mean(RootResp[Time_since_thaw == 0], na.rm = TRUE)) %>%
  filter(!is.na(month2_val), month2_val > 0) %>%
  mutate(rel_change = RootResp / month2_val) %>%
  ungroup()

rel_decline_root_derived %>% group_by(Locality, TempFactor) %>% summarize(month2_val = first(month2_val)) %>% arrange(month2_val)

rel_decline_root_derived %>% group_by(Time_since_thaw) %>%
  summarize(mean_rel = mean(rel_change, na.rm=TRUE), se=sd(rel_change,na.rm=TRUE)/sqrt(n()), n=n(), .groups="drop")

rel_decline_root_derived_clean <- rel_decline_root_derived %>%
  transmute(Locality, Time_since_thaw, SoilResp = RootResp, rel_change, TempFactor)

rel_decline_root_stable <- bind_rows(
  TrRoots_cc %>% filter(Temp==-5) %>% transmute(Locality, Time_since_thaw=round((Time-60.88)/30.44), SoilResp=RootResp, TempFactor="stable_-5"),
  TrRoots_cc %>% filter(Temp==5)  %>% transmute(Locality, Time_since_thaw=round((Time-60.88)/30.44), SoilResp=RootResp, TempFactor="stable_5")
) %>%
  group_by(TempFactor, Locality) %>%
  mutate(month2_val = mean(SoilResp[Time_since_thaw==0], na.rm=TRUE)) %>%
  filter(!is.na(month2_val), month2_val>0) %>%
  mutate(rel_change = SoilResp/month2_val) %>%
  ungroup()

rel_all_root <- bind_rows(
  rel_decline_root_derived_clean,
  rel_decline_root_stable
)

rel_all_root %>% group_by(TempFactor, Time_since_thaw) %>%
  summarize(mean_rel = mean(rel_change, na.rm=TRUE), se=sd(rel_change,na.rm=TRUE)/sqrt(n()), .groups="drop")

rel_root_model <- lm(rel_change ~ TempFactor * Time_since_thaw, data = rel_all_root %>% filter(Time_since_thaw>0))
anova(rel_root_model)
summary(rel_root_model)

# ---- Another figure ----
pred_month2_FT_root <- bind_rows(pred_5_5_root_val, pred_18_5_root_val) |> 
  mutate(TempFactor = rep(c("-5/5", "-18/5"), each = 6)) 


RootRef <- TrRoots_cc |> group_by(Locality, Temp, MonthsInc) |> 
  summarise(yref = mean(RootResp, na.rm = T), yrefSD = sd(RootResp, na.rm = T)) |> 
  filter(Temp != 0) |> filter(MonthsInc != 3) |>  
  mutate(TempF = factor(Temp, levels = c(-5, 5), labels = c("-5˚C", '+5˚C')))

pred_month2_FT_root2 <- rbind(pred_month2_FT_root, pred_month2_FT_root)
pred_month2_FT_root2$MonthsInc <- rep(c(2, 4), each = nrow(pred_month2_FT_root)) 
pred_month2_FT_root2[(nrow(pred_month2_FT_root)+1):nrow(pred_month2_FT_root2), c(2:4)] <- NA

TrLongerFull |> filter(Origin == 'Roots' & StabTemp == 'No' & MonthsInc != 3) |> 
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
  theme(axis.title.x = element_blank(), legend.position = c(0.1, 0.9)) +
  scale_x_discrete(labels = c(expression(atop("2160 m", "scree")),
                              expression(atop("2160 m", "meadow")),
                              expression(atop("2130 m", "meadow")),
                              expression(atop("1780 m", "meadow")),
                              expression(atop("1470 m", "forest")),
                              expression(atop("1330 m", "forest"))))

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

# boot_coefs_N2 and ratio0_draws/rel_boot have different lengths -> resample to align
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


# ---- supplementary figure ----
trajDFrame_roots$Series <- c('-5/5')
trajDFrame_roots18$Series <- c('-18/5')

rbind(trajDFrame_roots, trajDFrame_roots18) |> group_by(Series, Scenario, Locality, Month) |> 
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
  theme(axis.title.x = element_blank()) +
  labs(y = c("Relative effect of alternating thermal regime (nondimensional)"),
       color = c("Time (months)"), fill = c("Time (months)"))


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
