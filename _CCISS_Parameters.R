
# Setting applied to all runs

GCMs <- c("ACCESS1-0", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", 
          "GFDL-CM3", "GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", 
          "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")

edatopes <- c("B2", "C4", "D6")

spps.lookup <- read.csv("inputs/Tree speciesand codes_2.0_2May2019.csv")

edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")

BGCcolors.BC <- read.csv("inputs/BGCzone_Colorscheme.csv")

BGCcolors <- read.csv("inputs/WNAv11_Zone_Colours.csv")

BGCcolors.subzone <- read.csv("inputs/WNAv11_Subzone_Colours.csv")

rcps <- c("rcp45", "rcp85")

proj.years <- c(2025, 2055, 2085)

hist.years <- c(1995, 2004, 2005, 2009, 2014, 2018)

hist.year.name <- c("1991-2000", "1991-2017", "2001-2010", "2001-2018","2011-2018", "2018")