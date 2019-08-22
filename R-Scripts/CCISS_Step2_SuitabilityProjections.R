##======================================================================================
## CCISS Publication Scripts
## Step 2 - Species suitability projections
##======================================================================================

# Colin Mahony
# c_mahony@alumni.ubc.ca
# 778-288-4008
# July 21, 2019

# What this script does: calculate suitability outputs (but not stocking standards) equivalent to the CCISS tool, but simplified for fast running. 
# reason for script: the full CCISS script takes ~10-20 seconds per POI, meaning weeks of running time for a publishable raster grid. 
# Objective: produce results for a ~250,000-cell grid at a run time of maximum 8 hours.  
# general approach: this script achieves speed by: 
  # 1. running only for selected edatopes
  # 2. running on the whole grid vector, not looping on the POIs. 
  # 3. outputs and inputs in many small tables, rather than a few big tables. 
  # 4. running only for selected species

# Analysis Steps: 
#   /Edatope/RCP/Spp/proj.year/GCM: 
#     1. Project BGC unit 
#     2. find the site series associated with the edatope
#     3. obtain suitability of spp in site series

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

setwd("C:\\Colin\\Projects\\2019_CCISS")

grid <- "BC2kmGrid"

GCMs <-  c("ACCESS1-0","CanESM2","CCSM4","CESM1-CAM5","CNRM-CM5","CSIRO-Mk3-6-0", "GFDL-CM3","GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", "MIROC5", "MPI-ESM-LR","MRI-CGCM3")
rcps <- c("rcp45", "rcp85")
proj.years <- c(2025, 2055, 2085)
hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)
edatopes<- c("B2", "C4", "D6")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")
BGCcolors <- read.csv("C:\\Colin\\Projects\\2019_CCISS\\InputData\\BGCzone_Colorscheme.csv")

# Knowledge Tables
treesuit="TreeSpp_ESuit_v11_18"
SiteSeries_Use <-read.csv(paste("InputData/","SiteSeries_Use_5",".csv",sep=""),stringsAsFactors=FALSE,na.strings=".")
spps.lookup <- read.csv("InputData\\Tree speciesand codes_2.0_2May2019.csv")


#===============================================================================
# generate the vector of mapped BGCs
#===============================================================================

points <- read.csv(paste("InputData\\",grid,".csv", sep=""))
BGC <- points$ID2
BGC <- gsub(" ","",BGC)  
sort(table(BGC))

#===============================================================================
# Import BGC projections for each period
#===============================================================================

## reference period BGC
BGC.pred.ref <- as.character(read.csv(paste("OutputData\\BGC.pred", grid, "ref.csv", sep="."))[,1])

# Historical BGC
for(hist.year in hist.years){
  BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred", grid,hist.year,"csv", sep="."))[,1])
  assign(paste("BGC.pred",hist.year,sep="."), BGC.pred)
  print(hist.year)
}

# Future BGC
PredSum <- data.frame()
for(rcp in rcps){
  for(proj.year in proj.years){
    for(GCM in GCMs){
      BGC.pred <- as.character(read.csv(paste("OutputData\\BGC.pred",grid, GCM, rcp, proj.year,".csv", sep=""))[,1])
      assign(paste("BGC.pred", GCM, rcp, proj.year, sep="."), BGC.pred)
      # print(GCM)
      PredSum <- rbind(PredSum, as.data.frame(table(BGC.pred)))
    }
    print(proj.year)
  }
  print(rcp)
}
PredSum <- aggregate(PredSum$Freq, by=list(PredSum$BGC.pred), FUN=sum)
names(PredSum) <- c("BGC", "count")

# #===============================================================================
# # make a lookup table for the site series associated with the selected edatope in each BGC unit
# # the idea is to focus on a representative set of edatopes: subxeric [2B], mesic [4C] and hygric [6D]
# # NOTE: THIS IS RELIC CODE REPLACED BY IMPORTING WILL'S SITE SERIES LOOKUP TABLE BELOW
# #===============================================================================
#
# # Import edatopic table
# edatopename="Edatopic_v11_6"
# E1 <-read.csv(paste("InputData/",edatopename,".csv",sep=""),stringsAsFactors=FALSE,na.strings=".")
# E1 <- unique(E1)
# table(E1[,1])
# temp <- unique(E1[,1:2])
# length(temp[,1])
# length(unique(temp))
# #function for selecting n characters from the right of a string
# substrRight <- function(x, n){
#   substr(x, nchar(x)-n+1, nchar(x))
# }
# 
# edatopes<- c("B2", "C4", "D6")
# SiteLookup <- as.data.frame(matrix(rep(NA, length(unique(E1$MergedBGC))*(length(edatopes)+1)),length(unique(E1$MergedBGC)),length(edatopes)+1))
# names(SiteLookup) <- c("MergedBGC", edatopes)
# SiteLookup$MergedBGC <- unique(E1$MergedBGC)
# for(edatope in edatopes){
#   #subset the edatopic lookup for just the edatope of interest
#   temp <- E1[which(E1$Edatopic==edatope),]
#   temp <- temp[which(temp$Codes==""),-c(5:6)] #remove special sites
#   zonal <- which(substrRight(temp$SS_NoSpace,2)=="01")
#   if(edatope=="C4"){ # if zonal, populate zonal sites first
#     SiteLookup[,which(edatopes==edatope)+1] <- temp$SS_NoSpace[zonal][match(SiteLookup$MergedBGC, temp$MergedBGC[zonal])]
#     # then fill in with non-zonal where needed.
#     if(length(zonal)>0) temp <- temp[-zonal,]
#     SiteLookup[is.na(SiteLookup[,which(edatopes==edatope)+1]),which(edatopes==edatope)+1] <- temp$SS_NoSpace[match(SiteLookup$MergedBGC[is.na(SiteLookup[,which(edatopes==edatope)+1])], temp$MergedBGC)]
#   } else { ## if non-zonal, populate with non-zonal sites first
#     SiteLookup[,which(edatopes==edatope)+1] <- temp$SS_NoSpace[-zonal][match(SiteLookup$MergedBGC, temp$MergedBGC[-zonal])]
#     # then fill in with zonal where needed.
#     if(length(zonal)>0) temp <- temp[zonal,]
#     SiteLookup[is.na(SiteLookup[,which(edatopes==edatope)+1]),which(edatopes==edatope)+1] <- temp$SS_NoSpace[match(SiteLookup$MergedBGC[is.na(SiteLookup[,which(edatopes==edatope)+1])], temp$MergedBGC)]
#   }
#   ## then fill in with special sites where needed.
#   temp <- E1[which(E1$Edatopic==edatope),]
#   special <- -which(temp$Codes=="")
#   if(length(special)>0) temp <- temp[special,]
#   SiteLookup[is.na(SiteLookup[,which(edatopes==edatope)+1]),which(edatopes==edatope)+1] <- temp$SS_NoSpace[match(SiteLookup$MergedBGC[is.na(SiteLookup[,which(edatopes==edatope)+1])], temp$MergedBGC)]
# }
# SiteLookup
# SiteLookup[, which(names(SiteLookup)==edatope)] # XXX Note that there are a few NAs in here. need to get these out.
# dim(SiteLookup)

#===============================================================================
# import will's lookup table for the site series associated with the selected edatope in each BGC unit
#===============================================================================

SiteLookup <- data.frame(MergedBGC=unique(SiteSeries_Use$MergedBGC))
edatopes<- c("B2", "C4", "D6")
for(edatope in edatopes){
  # SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[match(SiteLookup[,1], SiteSeries_Use$MergedBGC[which(SiteSeries_Use$Use==edatope)])])
  SiteLookup <- cbind(SiteLookup, SiteSeries_Use$SS_NoSpace[which(SiteSeries_Use$Edatopic==edatope)])
  names(SiteLookup)[which(edatopes==edatope)+1] <- edatope
  }
str(SiteLookup)

#===============================================================================
# find the species suitability each projection/edatope/species combination
#===============================================================================

hist.years <- c(1995, 2004, 2005, 2009, 2014, 2017)

# Import suitability tables
S1 <- read.csv(paste("InputData/",treesuit,".csv",sep=""),stringsAsFactors=F,na.strings=".")
S1 <- unique(S1)

## EDA: are there suitabilities for all projected units? 
NoSuit <- PredSum[-which(PredSum$BGC%in%S1$BGC),]
NoSuit[rev(order(NoSuit$count)),]

## EDA: Which site series are missing suitabilities? 
for(edatope in edatopes) assign(paste("NoSuit", edatope, sep="."), SiteLookup[-which(SiteLookup[,which(names(SiteLookup)==edatope)]%in%S1$Unit),which(names(SiteLookup)==edatope)])
for(edatope in edatopes) print(get(paste("NoSuit", edatope, sep=".")))

## EDA: are there any units missing from the SiteSeries_Use table? 
BGClist <- unique(S1$BGC)
BGClist[-which(BGClist%in%SiteLookup$MergedBGC)]

# select the species to run the analysis on
spps <- unique(S1$Spp)
spps <- spps[-which(spps=="X")]
spps.candidate <- spps.lookup$TreeCode[-which(spps.lookup$Exclude=="x")]
spps <- spps[which(spps%in%spps.candidate)] 

for(spp in spps){
  for(edatope in edatopes){
    # get the suitability for the reference period predicted BGC.
    BGC.pred <- as.character(BGC.pred.ref) # get the BGC prediction

    # get the suitability for the selected species associated with each site series
    suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
    Suit.ref <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
    Suit.ref[is.na(Suit.ref)] <- 5 #set the NA values to suitability 5
    Suit.ref[Suit.ref==4] <- 5 #set 4 to suitability 5
    write.csv(Suit.ref, paste("OutputData\\Suit.ref", grid, spp, edatope, "csv", sep="."), row.names = F)

    for(hist.year in hist.years){
      BGC.pred <- as.character(get(paste("BGC.pred", hist.year, sep=".")))
      ## identify cells with no suitability interpretations
      bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
      bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]
      
      suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
      temp <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
      temp[is.na(temp)] <- 5 #set the NA values to suitability 5
      temp[temp==4] <- 5 #set 4 to suitability 5
      temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
      assign(paste("Suit", hist.year, sep="."), temp)
      write.csv(temp, paste("OutputData\\Suit", grid, hist.year, spp, edatope, "csv", sep="."), row.names = F)
    # print(hist.year)
      }
    
    # get the suitability for future periods, for each projection/edatope/species combination
    for(GCM in GCMs){
      for(rcp in rcps){
        for(proj.year in proj.years){
          # get the BGC projection and sub in the crosswalk between the modeled units and the table units
          BGC.pred <- as.character(get(paste("BGC.pred", GCM, rcp, proj.year, sep=".")))
          # BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)] <- as.character(Crosswalk$Tables[match(BGC.pred[which(BGC.pred%in%Crosswalk$Modeled)], Crosswalk$Modeled)]) # XXX THIS IS NOT CORRECT. NEED TO FIGURE OUT HOW TO INCORPORATE THE CROSSWALK TABLE PROPERLY. sub in the crosswalk between the modeled units and the table units

          ## identify cells with no suitability interpretations
          bgc.exotic <- (1:length(BGC.pred))[-which(BGC.pred%in%unique(BGC))]
          bgc.exotic.noSuit <- bgc.exotic[-which(BGC.pred[bgc.exotic]%in%unique(S1$BGC))]

          # get the suitability for the selected species associated with each site series
          suit <- S1$ESuit[which(S1$Spp==spp)][match(as.vector(unlist(SiteLookup[which(names(SiteLookup)==edatope)])), S1$Unit[which(S1$Spp==spp)])]
          temp <- suit[match(BGC.pred, SiteLookup$MergedBGC)]
          temp[is.na(temp)] <- 5 #set the NA values to suitability 5 (weights unsuitable a bit more heavily than suitable classes during averaging)
          temp[temp==4] <- 5 #set 4 to suitability 5
          temp[bgc.exotic.noSuit] <- NA # set cells with no suitabilty interpretatoin to NA
          assign(paste("Suit", GCM, rcp, proj.year, sep="."), temp)
          write.csv(temp, paste("OutputData\\Suit", grid, GCM, rcp, proj.year, spp, edatope, "csv", sep="."), row.names = F)
          # print(proj.year)
        }
        # print(rcp)
      }
      # print(GCM)
    }
    # print(edatope)
  }
  print(paste(spp, " (", round(which(spps==spp)/length(spps)*100, 0), "%)", sep=""))
}
