
# Setting applied to all runs

GCMs <- c("ACCESS1-0", "CanESM2", "CCSM4", "CESM1-CAM5", "CNRM-CM5", "CSIRO-Mk3-6-0", 
          "GFDL-CM3", "GISS-E2R", "HadGEM2-ES", "INM-CM4", "IPSL-CM5A-MR", "MIROC-ESM", 
          "MIROC5", "MPI-ESM-LR", "MRI-CGCM3")

Columns <- unique(c("PPT05", "PPT06", "PPT07", "PPT08", "PPT09", "PPT_at", 
                    "PPT_wt", "CMD07", "CMD", "MAT", "PPT_sm", "Tmin_wt", "Tmax_sm",
                    vars[!vars %in% c("PPT_MJ", "PPT_JAS", "PPT.dormant", "CMD.def", "CMDMax", "CMD.total")]))

rcps <- c("rcp45", "rcp85")

proj.years <- c(2025, 2055, 2085)

edatopes <- c("B2", "C4", "D6")

spps.lookup <- fread("inputs/Tree speciesand codes_2.0_2May2019.csv")

edatope.name <- c("Subxeric-poor", "Mesic-medium", "Hygric-rich")

BGCcolors <- fread("inputs/BGCzone_Colorscheme.csv")
