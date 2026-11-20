#===============================================================================
#                               Tatra Mountains 
#                             Incubation experiment
#
# This script loads original data from xlsx files and produce reproducible 
# source csv file  
################################################################################

#===============================================================================
#Libraries and ggplot theme are loaded from separate files
#===============================================================================
mytheme <- readRDS('ggtheme.rds')
source("packages.R")
#===============================================================================
# Data
#===============================================================================
Tr <- read.xlsx('SourceDataXLSXL/Tatry_valecky_teplotni-inkubace2024_GC.xlsx', 8)
TrPudni <- read.xlsx('SourceDataXLSX/valecky_teplota.xlsx', 1)
################################################################################


#===============================================================================
# Root and soil-derived respiration - two-pool mixing model
#===============================================================================
# ---- Root-free (control) and root amended data sets are filtered ----
TrControl <- subset(Tr, Roots == 'No')
TrRoots <- subset(Tr, Roots == 'Yes')

# Matching the order to make sure that calculation is correct
TrControl <- TrControl |> arrange(MonthsInc, Locality, TempFactor, d13Ccalc)
TrRoots <- TrRoots |> arrange(MonthsInc, Locality, TempFactor, d13Ccalc)

TrRoots$d13Ccontrol <- TrControl$d13Ccalc
TrRoots$RespControl <- TrControl$Resp

# Missing data are replaced by the mean of other values
for(i in unique(Tr$MonthsInc)){
  for(n in unique(Tr$Locality)){
    for(l in unique(Tr$TempFactor)){
      if(sum(is.na(TrRoots[(TrRoots$MonthsInc == i & 
                            TrRoots$Locality == n &
                            TrRoots$TempFactor == l), 'd13Ccalc'])) > 0 &
         sum(is.na(TrRoots[(TrRoots$MonthsInc == i & 
                            TrRoots$Locality == n &
                            TrRoots$TempFactor == l), 'd13Ccontrol'])) == 0){
        TrRoots[(TrRoots$MonthsInc == i & 
                   TrRoots$Locality == n &
                   TrRoots$TempFactor == l &
                   is.na(TrRoots$d13Ccalc)), 'd13Ccalc'] <- as.numeric(mean(TrRoots[(TrRoots$MonthsInc == i & 
                                                                            TrRoots$Locality == n &
                                                                            TrRoots$TempFactor == l), 'd13Ccalc'], na.rm = T))
      } else {
        if(sum(is.na(TrRoots[(TrRoots$MonthsInc == i & 
                              TrRoots$Locality == n &
                              TrRoots$TempFactor == l), 'd13Ccontrol'])) > 0 &
           sum(is.na(TrRoots[(TrRoots$MonthsInc == i & 
                              TrRoots$Locality == n &
                              TrRoots$TempFactor == l), 'd13Ccalc'])) == 0){
          TrRoots[(TrRoots$MonthsInc == i & 
                     TrRoots$Locality == n &
                     TrRoots$TempFactor == l &
                     is.na(TrRoots$d13Ccontrol)), 'd13Ccontrol'] <- as.numeric(mean(TrRoots[(TrRoots$MonthsInc == i & 
                                                                              TrRoots$Locality == n &
                                                                              TrRoots$TempFactor == l), 'd13Ccontrol'], na.rm = T))
        }
      }
    }
  }
}


# Order again  
TrRoots <- TrRoots |> arrange(MonthsInc, Locality, TempFactor, d13Ccalc)

# ---- Fraction of root-derived respiration ----
TrRoots$frRoot <- with(TrRoots, (d13Ccalc - d13Ccontrol)/(d13CRoots - d13Ccontrol))
summary(TrRoots$frRoot)
sum((TrRoots$frRoot < 0), na.rm = T)
sum((TrRoots$frRoot > 1), na.rm = T)

TrRoots$frRoot <- ifelse(TrRoots$frRoot<0, 0, TrRoots$frRoot)
TrRoots$frRoot <- ifelse(TrRoots$frRoot>1, 1, TrRoots$frRoot)

# ---- Pasting to original data frame Tr ----
identical(TrRoots$Locality, Tr[Tr$Roots == "No", 'Locality'])

TrRoots <- TrRoots |> arrange(Locality, TempFactor, MonthsInc)
Tr <- Tr |> arrange(Roots, Locality, TempFactor, MonthsInc)

identical(TrRoots$Locality, Tr[Tr$Roots == "Yes", 'Locality'])

Tr$frRoot <- ifelse(Tr$Roots == 'No', 0, TrRoots$frRoot)
Tr$SoilResp <- with(Tr, Resp*(1 - frRoot))
Tr$RootResp <- with(Tr, Resp*frRoot)

# ---- Reformatting data frame to long format ---- 
# a new factor 'origin' with levels "Control soil", "Soil with roots" and "Roots" is created - this will be handy for data subsetting
TrLonger <- rbind(Tr, Tr[271:540, ])
TrLonger$RespAll <- c(Tr$SoilResp, Tr$RootResp[271:540])
TrLonger$Origin <- rep(c("Control soil", 'Soil with roots', 'Roots'), each = 270)

# ---- Merging with soil data ----
colnames(TrPudni)[1:2] <- c('Locality', 'Roots')
colnames(TrPudni)
unique(Tr$TempFactor)
TrPudni$TempFactor <- character(length = nrow(TrPudni))
for(i in 1:nrow(TrPudni)){
  if(TrPudni$temperature[i] == 'A'){
    TrPudni$TempFactor[i] <- '5'
  }else{
    if(TrPudni$temperature[i] == 'B'){
      TrPudni$TempFactor[i] <- '0'
    }else{
      if(TrPudni$temperature[i] == 'C'){
        TrPudni$TempFactor[i] <- '-5'
      }else{
        if(TrPudni$temperature[i] == 'D'){
          TrPudni$TempFactor[i] <- '-5/5'
        }else{
          TrPudni$TempFactor[i] <- '-18/5'
        }
      }
    }
  }
}
TrPudni$Locality <- as.factor(TrPudni$Locality)
TrPudniLonger <- rbind(TrPudni, TrPudni[c(91:180), ])
TrPudniLonger$Origin <- rep(unique(TrLonger$Origin), each = 90)
TrPudniLonger$triplicate <- rep(c(1,2,3), times = 90)
TrPudniLonger$key <- paste0(TrPudniLonger$Locality, '-', TrPudniLonger$TempFactor, '-',
                            TrPudniLonger$Origin, '-', TrPudniLonger$triplicate)
sum(duplicated(TrPudniLonger$key)) 
TrLonger$triplicate <- rep(c(1,2,3), times = 90*3)
TrLonger$key <- paste0(TrLonger$Locality, '-', TrLonger$TempFactor, '-',
                       TrLonger$Origin, '-', TrLonger$triplicate)
sum(duplicated(TrLonger$key)) 

TrLongerFull <- merge(TrLonger, TrPudniLonger, by = 'key', no.dups = F)

# Initial root C concentration 
TrLongerFull$C0 <- ifelse(TrLongerFull$Origin.x == 'Roots', 100*0.55/40/TrLongerFull$susina, TrLongerFull$soil_C)

# Removing duplicates
names(TrLongerFull)

TrLongerFull <- TrLongerFull[, -c(1, 20:22, 55:57)]
names(TrLongerFull)
colnames(TrLongerFull)[c(1, 5, 6, 18)] <- c("Locality", "TempFactor", "Roots", "Origin") 

#===========================Exporting - this is the working data set
write.csv(TrLongerFull, 'TrLongerFull.csv', row.names = F)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF SCRIPT
