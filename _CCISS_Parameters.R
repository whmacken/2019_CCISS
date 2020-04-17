
# Setting applied to all runs

grid <- "BC2kmGrid"

GCMs <- c("ACCESS1-0", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", 
          "GFDL-CM3", "GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", 
          "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")

edatopes <- c("B2", "C4", "D6")
edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")

rcps <- c("rcp45", "rcp85")
rcp.name=c("RCP4.5", "RCP8.5")

proj.years <- c(2025, 2055, 2085)
proj.year.name=c("2020s", "2050s", "2080s")

hist.years <- c(1995, 2004, 2005, 2009, 2014, 2018)
hist.year.name <- c("1991-2000", "1991-2018", "2001-2010", "2001-2018","2011-2018", "2018")

model <- "35Var_6_2"
fname <- "inputs/models/WNAv11_35_VAR_SubZone_ranger.Rdata"

####Lookup tables

BGCcolors.BC <- read.csv("lookup/BGCzone_Colorscheme.csv")
BGCcolors <- read.csv("lookup/WNAv11_Zone_Colours.csv")
BGCcolors.subzone <- read.csv("lookup/WNAv11_Subzone_Colours.csv")
treesuit <- read.csv("lookup/ESuitv11_21.csv")
SiteSeries_Use <-read.csv("lookup/SiteSeries_Use_CCISSpaper_13Apr2020.csv",stringsAsFactors=FALSE,na.strings=".")
spps.lookup <- read.csv("lookup/Tree speciesand codes_2.0_2May2019.csv")
BGCs_notin_THLB <- read.csv("lookup/BGCs_notin_THLB.csv")

#BGC zone color scheme
BGCcolors$colour <- as.character(BGCcolors$colour)
BGCcolors$colour[match(BGCcolors.BC$zone, BGCcolors$classification)] <- as.character(BGCcolors.BC$HEX)
ColScheme.zone <- factor(BGCcolors$colour, levels=BGCcolors$colour)
zones <- factor(BGCcolors$classification, levels=BGCcolors$classification)
