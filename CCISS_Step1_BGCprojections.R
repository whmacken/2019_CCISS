
##======================================================================================
## CCISS Publication Scripts
## Step 1 - BGC projections for BC study area
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019

require (RGtk2)
require(plyr)
require (rChoiceDialogs)
require (data.table)
require(doBy)
require (utils)
require(labdsv)
require(tools )
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


# rm(list=ls())

#===============================================================================
# Set analysis Parameters
#===============================================================================

setwd("C:/GitHub/2019_CCISS/")

grid <- "BC2kmGrid"

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")    

###Load random forest model
model = "5.1"
fname="InputData\\BGCv11_AB_USA_16VAR_SubZone_RFmodel.Rdata"
load(fname)
# rownames(importance(BGCmodel))


GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
edatopes<- c("B2", "C4", "D6")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0_2May2019.csv")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")


#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- points$ID2
table(BGC)

BGC <- gsub(" ","",BGC)  
table(BGC)
sort(table(BGC))

# # merge small units
# BGC[which(BGC=="CMAwh")] <- "CMAun"
# BGC[which(BGC=="MHun")] <- "MHunp"
# BGC[which(BGC=="ESSFxvw")] <- "ESSFxvp"
# BGC[which(BGC=="ESSFdcp")] <- "ESSFdcw"

#BGC zones
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")
zone <- rep(NA, length(BGC))
for(i in BGCcolors$zone){ zone[grep(i,BGC)] <- i }
table(zone)


# #===============================================================================
# # select variables from ClimateWNA file and subset by GCM
# # NB: is commented out because you only need to do this once: it just reduces ram needed to load in the ClimateBC data
# #===============================================================================
# 
# setwd("C:\\Colin\\Projects\\2019_CCISS")
# Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", "PPT_wt", "CMD07", "CMD", "MAT", "PPT_sm", "Tmin_wt", "Tmax_sm",  rownames(importance(BGCmodel))[-which(rownames(importance(BGCmodel))%in%c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total"))]))
# 
# #first batch of 8 models
# fplot=paste("InputData\\", grid, "_48GCMs_MSYT.csv", sep="")
# Y1 <- fread(fplot, select = c("GCM", Columns), stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
# models <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R")
# for(model in models){
#   temp <- Y1[grep(model,Y1$GCM),]
#   write.csv(temp, paste("InputData\\", grid, "_", model, "_BioVars.csv", sep=""), row.names = F)
#   print(which(models==model))
# }
# 
# #second batch of 7 models
# fplot=paste("InputData\\", grid, "_42GCMs_MSYT.csv", sep="")
# Y1 <- fread(fplot, select = c("GCM", Columns), stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
# models <-  c("HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
# for(model in models){
#   temp <- Y1[grep(model,Y1$GCM),]
#   write.csv(temp, paste("InputData\\", grid, "_", model, "_BioVars.csv", sep=""), row.names = F)
#   print(which(models==model))
# }


#===============================================================================
# BGC Projections for reference period
#===============================================================================

setwd("C:/GitHub/2019_CCISS")
Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", "PPT_wt", "CMD07", "CMD", "MAT", "PPT_sm", "Tmin_wt", "Tmax_sm",  rownames(importance(BGCmodel))[-which(rownames(importance(BGCmodel))%in%c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total"))]))

fplot=paste("InputData\\", grid, "_Normal_1961_1990MSY.csv", sep="")

Y0 <- fread(fplot, select=Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv

Y0 <- Y0[!is.na(Y0[,2]),]

#####generate some additional variables
Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
Y0$CMDMax <- Y0$CMD07
Y0$CMD.total <- Y0$CMD.def + Y0$CMD

##Predict future subzones######
BGC.pred.ref <- predict(BGCmodel, Y0)
write.csv(BGC.pred.ref, paste("OutputData\\BGC.pred.ref", grid, "csv", sep="."), row.names = F)

## Write Climate file ######
write.csv(Y0, paste("InputData\\", grid, "_1961_1990_BioVars.csv", sep=""), row.names = F)

#===============================================================================
# BGC Projections for historical decades
#===============================================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")
hist.years <- c(1995, 2005)
hist.periods <- c("1991_2000", "2001_2010")

for(hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years==hist.year)]
  fplot=paste("InputData\\", grid, "_Decade_", hist.period, "MSY.csv", sep="")

  Y0 <- fread(fplot, select=Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
  Y0 <- Y0[!is.na(Y0[,2]),]
# str(Y0)

  #####generate some additional variables
  Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
  Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
  Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
  Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
  Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
  Y0$CMDMax <- Y0$CMD07
  Y0$CMD.total <- Y0$CMD.def + Y0$CMD

  ##Predict future subzones######
  assign(paste("BGC.pred",hist.year,sep="."), predict(BGCmodel, Y0))
  write.csv(get(paste("BGC.pred",hist.year,sep=".")), paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."), row.names = F)

  ## Write Climate file ######
  write.csv(Y0, paste("InputData\\", grid, "_",hist.year,"_BioVars.csv", sep=""), row.names = F)

  print(hist.year)
}

#===============================================================================
# BGC Projections for last decade
#===============================================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")
  fplot=paste("InputData\\", grid, "_2011-2017MSYT.csv", sep="")

  Y0 <- fread(fplot, select=c("ID1", "Year", Columns), stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
  # Y0 <- Y0[!is.na(Y0[,2]),]
str(Y0)

  #####generate some additional variables
  Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
  Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
  Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
  Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
  Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
  Y0$CMDMax <- Y0$CMD07
  Y0$CMD.total <- Y0$CMD.def + Y0$CMD

  # Extract the year 2017
  Y2017 <- Y0[which(Y0$Year==2017),]
str(Y2017)

  # Calculate mean of 2011-2017 period
  Y1 <- aggregate(Y0, by=list(Y0$ID1), FUN=mean, na.rm=T)[,-1]
  Y1 <- Y1[match(Y2017$ID1, Y1$ID1),]

  ##Predict BGC units######
  BGC.pred.2017 <- predict(BGCmodel, Y2017)
  BGC.pred.2014 <- predict(BGCmodel, Y1)
  write.csv(BGC.pred.2017, paste("OutputData\\BGC.pred", grid,"2017.csv", sep="."), row.names = F)
  write.csv(BGC.pred.2014, paste("OutputData\\BGC.pred", grid,"2014.csv", sep="."), row.names = F)

  ## Write Climate file ######
  write.csv(Y2017, paste("InputData\\", grid, "_2017_BioVars.csv", sep=""), row.names = F)
  write.csv(Y1, paste("InputData\\", grid, "_2014_BioVars.csv", sep=""), row.names = F)

#===============================================================================
# BGC Projections for other historical normals
#===============================================================================

setwd("C:\\Colin\\Projects\\2019_CCISS")
hist.years <- c(1995, 2005)
hist.periods <- c("1991_2000", "2001_2010")

mean(c(1991,2017))
mean(c(2001,2017))

# read in the data for the 1990s and 2000s
for(hist.year in hist.years){
  hist.period <- hist.periods[which(hist.years==hist.year)]
  fplot=paste("InputData\\", grid, "_Decade_", hist.period, "MSY.csv", sep="")

  Y0 <- fread(fplot, select=Columns, stringsAsFactors = FALSE, data.table = FALSE) #fread is faster than read.csv
  Y0 <- Y0[!is.na(Y0[,2]),]

  #####generate some additional variables
  Y0$PPT_MJ <- Y0$PPT05 + Y0$PPT06 # MaY/June precip
  Y0$PPT_JAS <- Y0$PPT07 + Y0$PPT08 + Y0$PPT09 # July/Aug/Sept precip
  Y0$PPT.dormant <- Y0$PPT_at + Y0$PPT_wt # for calculating spring deficit
  Y0$CMD.def <- 500 - (Y0$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
  Y0$CMD.def [Y0$CMD.def < 0] <- 0 #negative values set to zero = no deficit
  Y0$CMDMax <- Y0$CMD07
  Y0$CMD.total <- Y0$CMD.def + Y0$CMD
  assign(paste("Y",hist.year,sep="."), Y0)
  print(hist.year)
}

# 2011-2017 period already exported in previous phase, so read that in
Y.2014 <- read.csv(paste("InputData\\", grid, "_2014_BioVars.csv", sep=""))
Y.2014 <- Y.2014[,-which(names(Y.2014)%in%c("ID1", "Year"))]
str(Y.2014)

Y.2009 <- (Y.2005*10+Y.2014*7)/17
Y.2004 <- (Y.1995*10+Y.2005*10+Y.2014*7)/27
str(Y.2009)
str(Y.2004)

hist.years <- c(2004, 2009)
for(hist.year in hist.years){

  ##Predict future subzones######
  assign(paste("BGC.pred",hist.year,sep="."), predict(BGCmodel, get(paste("Y",hist.year,sep="."))))
  write.csv(get(paste("BGC.pred",hist.year,sep=".")), paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."), row.names = F)

  ## Write Climate file ######
  write.csv(get(paste("Y",hist.year,sep=".")), paste("InputData\\", grid, "_",hist.year,"_BioVars.csv", sep=""), row.names = F)

  print(hist.year)
}

#===============================================================================
# BGC Projections for future periods
#===============================================================================
rcp="rcp45"

for(GCM in GCMs){
  Y1 <- fread(paste("InputData\\", grid, "_", GCM, "_BioVars.csv", sep=""), stringsAsFactors = FALSE, data.table = FALSE)
  # Y1 <- Y1[!is.na(Y1[,2]),]

  #####generate some additional variables
  Y1$PPT_MJ <- Y1$PPT05 + Y1$PPT06 # MaY/June precip
  Y1$PPT_JAS <- Y1$PPT07 + Y1$PPT08 + Y1$PPT09 # July/Aug/Sept precip
  Y1$PPT.dormant <- Y1$PPT_at + Y1$PPT_wt # for calculating spring deficit
  Y1$CMD.def <- 500 - (Y1$PPT.dormant)# start of growing season deficit original value was 400 but 500 seems better
  Y1$CMD.def [Y1$CMD.def < 0] <- 0 #negative values set to zero = no deficit
  Y1$CMDMax <- Y1$CMD07
  Y1$CMD.total <- Y1$CMD.def + Y1$CMD

  ## assign single vectors to RCPs and proj.years
  Ystr <- strsplit(Y1[,1], "_")
  Y4 <- matrix(unlist(Ystr), ncol=3, byrow=TRUE)

  proj.years <- unique(Y4[,3])
  rcps <- unique(Y4[,2])
  for(rcp in rcps){
    for(proj.year in proj.years){
      temp <- Y1[which(Y4[,2]==rcp & Y4[,3]==proj.year),]
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), predict(BGCmodel, temp))
      write.csv(get(paste("BGC.pred", GCM, rcp, proj.year, sep=".")), paste("OutputData\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""), row.names = F)
      print(proj.year)
    }
    print(rcp)
  }
  print(GCM)
}


