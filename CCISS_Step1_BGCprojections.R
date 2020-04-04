
## ======================================================================================
## CCISS Publication Scripts Step 1 - BGC projections for BC study area
## ======================================================================================

# Colin Mahony c_mahony@alumni.ubc.ca 778-288-4008 July 21, 2019

require(RGtk2)
require(plyr)
require(rChoiceDialogs)
require(data.table)
require(doBy)
require(utils)
require(labdsv)
require(tools)
require(svDialogs)
require(tcltk)
require(randomForest)
require(foreach)
require(dplyr)
require(reshape2)
require(reshape)
library(doParallel)
require(data.table)
library(MASS)
library(scales)
library(stats)
library(rgl)
library(RColorBrewer)
library(FNN)
library(igraph)
library(raster)
library(maps)
library(mapdata)
library(maptools)
library(sp)
library(colorRamps)
library(rgeos)
library(rgdal)
library(foreign)
require(data.table)
require(ranger)
require(tidyverse)


rm(list = ls())

# ===============================================================================
# Set analysis Parameters
# ===============================================================================

grid <- "BC2kmGrid"

GCMs <- c("ACCESS1-0", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", 
          "GFDL-CM3", "GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", 
          "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")

addVars <- function(dat){
  dat$PPT_MJ <- dat$PPT05 + dat$PPT06  # MaY/June precip
  dat$PPT_JAS <- dat$PPT07 + dat$PPT08 + dat$PPT09  # July/Aug/Sept precip
  dat$PPT.dormant <- dat$PPT_at + dat$PPT_wt  # for calculating spring deficit
  dat$CMD.def <- 500 - (dat$PPT.dormant)  # start of growing season deficit original value was 400 but 500 seems better
  dat$CMD.def[dat$CMD.def < 0] <- 0  #negative values set to zero = no deficit
  dat$CMDMax <- dat$CMD07
  dat$CMD.total <- dat$CMD.def + dat$CMD
  X1$CMD.grow <- X1$CMD05 + X1$CMD06 +X1$CMD07 +X1$CMD08 +X1$CMD09
  X1$DD5.grow <- X1$DD5_05 + X1$DD5_06 + X1$DD5_07 + X1$DD5_08 + X1$DD5_09
  X1$CMDMax <- X1$CMD07 # add in so not removed below
  X1$DDgood <- X1$DD5 - X1$DD18
  X1$DDnew <- (X1$DD5_05 + X1$DD5_06 +X1$DD5_07  + X1$DD5_08)  - (X1$DD18_05 + X1$DD18_06 +X1$DD18_07 +X1$DD18_08)
  X1$TmaxJuly <- X1$Tmax07
  return(dat)
}


### Load random forest model
model <- "16Var_6_2"
fname <- "inputs/models/WNAv11_16_VAR_SubZone_ranger.Rdata"
load(fname)
#rownames(importance(BGCmodel)) ### shows the variable used in the RFmodel


GCMs <- c("ACCESS1-0", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", 
          "GFDL-CM3", "GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", 
          "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2018)
edatopes <- c("B2", "C4", "D6")
spps.lookup <- fread("inputs/Tree speciesand codes_2.0_2May2019.csv")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
BGCcolors <- fread("inputs/BGCzone_Colorscheme.csv")


# ===============================================================================
# BGC Projections for reference period
# ===============================================================================
vars <- as.data.frame(BGCmodel$variable.importance)
vars <- row.names(vars)
# setwd('C:/GitHub/2019_CCISS')
Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", 
                    "PPT_wt", "CMD07", "CMD", "MAT", "PPT_sm", "Tmin_wt", "Tmax_sm",
                    vars[!vars %in% c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total")]))

fplot <- paste("inputs/", grid, "_Normal_1961_1990MSY.csv", sep = "")

Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE)  #fread is faster than read.csv

Y0 <- Y0[!is.na(Y0[, 2]), ]

Y0 <- addVars(Y0)

Y0 <- Y0 %>% dplyr::select(vars)

## Predict future subzones######
BGC.pred.ref <- predict(BGCmodel, Y0)
#dir.create("./outputs")
fwrite(BGC.pred.ref$predictions, paste("outputs/BGC.pred.ref", grid, "csv", sep = "."), 
          row.names = F)

## Write Climate file ######
fwrite(Y0, paste("inputs/", grid, "_1961_1990_", model,".csv", sep = ""))

# ===============================================================================
# BGC Projections for historical decades
# ===============================================================================

# setwd('C:\\Colin\\Projects\\2019_CCISS')
hist.years <- c(1985, 1995, 2005, 2014)
hist.periods <- c("1981_1990", "1991_2000", "2001_2010", "2011_2018")
#hist.year="1995"

for (hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years == hist.year)]
  fplot <- paste("inputs/", grid, "_Decade_", hist.period, "MSY.csv", 
                 sep = "")
  
  Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE)  #fread is faster than read.csv
  Y0 <- Y0[!is.na(Y0[, 2]), ]
  # str(Y0)
  Y0 <- addVars(Y0)
  Y0 <- Y0 %>% dplyr::select(vars)
  
  ## Predict future subzones######
BGC.pred <- predict(BGCmodel, Y0)
  
#assign(paste("BGC.pred", hist.year, sep = "."), 
  fwrite(BGC.pred$predictions, paste0("./outputs/BGC.pred", grid, hist.year, ".csv", sep = "."), row.names = F)
  
  ## Write Climate file ######
  fwrite(Y0, paste("inputs/", grid, "_", hist.year, "_", model, ".csv",  sep = ""), row.names = F)
  
  print(hist.year)
}

# ===============================================================================
# BGC Projections for last year available
# ===============================================================================

fplot <- paste("inputs/", grid, "_2011-2018MSY.csv", sep = "")

Y0 <- fread(fplot, select = c("ID1", "Year", Columns), stringsAsFactors = FALSE, 
            data.table = FALSE)  #fread is faster than read.csv
# Y0 <- Y0[!is.na(Y0[,2]),]
str(Y0)

Y0 <- addVars(Y0)
Y0 <- Y0 %>% dplyr::select(vars)
# Extract the year 2018
Y0 <- Y0[which(Y0$Year == 2018), ]

## Predict BGC units######
BGC.pred.2018 <- predict(BGCmodel, Y0)
fwrite(BGC.pred.2018$predictions, paste("outputs/BGC.pred", grid, "2018.csv", sep = "."), 
          row.names = F)

## Write Climate file ######
fwrite(Y0, paste("inputs/", grid, "_2018_BioVars.csv", sep = ""))

# ===============================================================================
# BGC Projections for other historical normals
# ===============================================================================

# setwd('C:\\Colin\\Projects\\2019_CCISS')
hist.years <- c(1995, 2005)
hist.periods <- c("1991_2000", "2001_2010")

mean(c(1991, 2018))
mean(c(2001, 2018))

# read in the data for the 1990s and 2000s
for (hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years == hist.year)]
  fplot <- paste("inputs/", grid, "_Decade_", hist.period, "MSY.csv", 
                 sep = "")
  
  Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE)  #fread is faster than read.csv
  Y0 <- Y0[!is.na(Y0[, 2]), ]
  
  Y0 <- addVars(Y0)
  assign(paste("Y", hist.year, sep = "."), Y0)
  print(hist.year)
}

# 2011-2018 period already exported in previous phase, so read that in
Y.2014 <- fread(paste("inputs/", grid, "_2014_BioVars.csv", sep = ""))
Y.2014 <- Y.2014[, -which(names(Y.2014) %in% c("ID1", "Year"))]
str(Y.2014)

Y.2009 <- (Y.2005 * 10 + Y.2014 * 7)/17
Y.2004 <- (Y.1995 * 10 + Y.2005 * 10 + Y.2014 * 7)/27
str(Y.2009)
str(Y.2004)

hist.years <- c(2004, 2009)
for (hist.year in hist.years){
  
  ## Predict future subzones######
  BGC.pred <- predict(BGCmodel, get(paste("Y", hist.year, sep = ".")))
  fwrite(BGC.pred$predictions, paste("outputs/BGC.pred", grid, hist.year, "csv", sep = "."), row.names = F)
  
  ## Write Climate file ######
  fwrite(get(paste("Y", hist.year, sep = ".")), paste("inputs/", grid, 
                                                         "_", hist.year, "_BioVars.csv", sep = ""), row.names = F)
  
  print(hist.year)
}

# ===============================================================================
# BGC Projections for future periods
# ===============================================================================

library(tidyr)
Y0 <- fread(paste0("./inputs/",grid,"_90 GCMsMSY.csv", sep = ""), select = c("Year", "ID1","ID2", Columns))   ##
Y0 <- separate(Y0, Year, into = c("Model","Scenario","FuturePeriod"), sep = "_", remove = T)
Y0$FuturePeriod <- gsub(".gcm","",Y0$FuturePeriod)
Y0 <- Y0 %>% select(Year, ID1,ID2, vars)

require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-2)
cl <- makeCluster(coreNum)
registerDoParallel(cl, cores = coreNum)

out <- foreach(GCM = GCMs, .combine = rbind) %:%
  foreach(rcp = rcps, .combine = rbind) %:%
  foreach(proj.year = proj.years, .combine = rbind, .packages = c("data.table","ranger", "dplyr")) %do% {
    sub <- Y0[Y0$Model == GCM & Y0$Scenario == rcp & Y0$FuturePeriod == proj.year,]
    sub <- addVars(sub)
    subPred <- data.frame(ID = sub$ID1)
    subPred$BGC.pred <- predict(BGCmodel, sub)
    
    fwrite(subPred, paste("outputs/BGC.pred", grid, GCM, rcp, proj.year, ".csv", sep = ""))
    temp <- aggregate(ID ~ BGC.pred, data = subPred, FUN = length) %>% 
      mutate(Model = GCM, Scn = rcp, FuturePeriod = proj.year)
    temp
  }



