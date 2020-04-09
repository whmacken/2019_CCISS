
## ======================================================================================
## CCISS Publication Scripts Step 1 - BGC projections for BC study area
## ======================================================================================

# Colin Mahony c_mahony@alumni.ubc.ca 778-288-4008 July 21, 2019
rm(list = ls())
source("./_CCISS_Packages.R") ## packages required
source("./_CCISS_Functions.R") ## common functions
source("./_CCISS_Parameters.R") ## settings used through all scripts

# ===============================================================================
# Set analysis Parameters
# ===============================================================================

grid <- "BC2kmGrid"

### Load random forest model
model <- "18Var_6_2"
fname <- "inputs/models/WNAv11_18_VAR_SubZone_ranger.Rdata"
load(fname)
#rownames(importance(BGCmodel)) ### shows the variable used in the RFmodel
vars <- as.data.frame(BGCmodel$variable.importance)
vars <- row.names(vars)
# setwd('C:/GitHub/2019_CCISS')



# ===============================================================================
# BGC Projections for reference period
# ===============================================================================


fplot <- paste("inputs/", grid, "_Normal_1961_1990MSY.csv", sep = "")

Y0 <- fread(fplot, select = Columns, stringsAsFactors = FALSE, data.table = FALSE)  #fread is faster than read.csv

Y0 <- Y0[!is.na(Y0[, 2]), ]

Y0 <- addVars(Y0)

Y0 <- Y0 %>% dplyr::select(all_of(vars))

## Predict future subzones######
BGC.pred.ref <- predict(BGCmodel, Y0)
#dir.create("./outputs")
fwrite(list(BGC.pred.ref$predictions), paste("outputs/BGC.pred.ref", grid, "csv", sep = "."), 
          row.names = F)

## Write Climate file ######
# fwrite(Y0, paste("inputs/", grid, "_1961_1990_", model,".csv", sep = ""))

# ===============================================================================
# BGC Projections for historical decades
# ===============================================================================

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
  fwrite(list(BGC.pred$predictions), paste0("outputs/BGC.pred", grid, hist.year, ".csv", sep = "."), row.names = F)
  
  ## Write Climate file ######
  # fwrite(Y0, paste("inputs/", grid, "_", hist.year, "_", model, ".csv",  sep = ""), row.names = F)
  
  print(hist.year)
}

# ===============================================================================
# BGC Projections for last year available
# ===============================================================================

fplot <- paste("inputs/", grid, "_Year_2018MSY.csv", sep = "")

Y0 <- fread(fplot, select = c("ID1",  Columns), stringsAsFactors = FALSE, 
            data.table = FALSE)  #fread is faster than read.csv
# Y0 <- Y0[!is.na(Y0[,2]),]
str(Y0)

#Y0 <- Y0[which(Y0$Year == 2018), ]

Y0 <- addVars(Y0)
Y0 <- Y0 %>% dplyr::select(vars)
# Extract the year 2018

## Predict BGC units######
BGC.pred.2018 <- predict(BGCmodel, Y0)
fwrite(list(BGC.pred.2018$predictions), paste("outputs/BGC.pred", grid, "2018.csv", sep = "."), 
          row.names = F)

## Write Climate file ######
# fwrite(Y0, paste("inputs/", grid, "_2018_",model,".csv", sep = ""))

# ===============================================================================
# BGC Projections for other historical normals
# ===============================================================================

hist.years <- c(1995, 2005, 2014)
hist.periods <- c("1991_2000", "2001_2010", "2011_2018")

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

Y.2009 <- (Y.2005 * 10 + Y.2014 * 7)/18
Y.2004 <- (Y.1995 * 10 + Y.2005 * 10 + Y.2014 * 7)/28
str(Y.2009)
str(Y.2004)

hist.years <- c(2004, 2009)
for (hist.year in hist.years){
  
  ## Predict future subzones######
  BGC.pred <- predict(BGCmodel, get(paste("Y", hist.year, sep = ".")))
  fwrite(list(BGC.pred$predictions), paste("outputs/BGC.pred", grid, hist.year, "csv", sep = "."), row.names = F)
  
  ## Write Climate file ######
  # fwrite(get(paste("Y", hist.year, sep = ".")), paste("inputs/", grid, "_", hist.year, "_", model, ".csv", sep = ""), row.names = F)
  
  print(hist.year)
}

# ===============================================================================
# BGC Projections for future periods
# ===============================================================================

GCM=GCMs[1]
rcp=rcps[1]
proj.year=proj.years[1]
Columns <- c("Year",  Columns)
##===============
###Option if one big all GCM file
Y0 <- fread("./inputs/BC2kmGrid_90 GCMsMSY.csv", select = Columns, stringsAsFactors = FALSE, data.table = FALSE)
Y0 <- separate(Y0, Year, into = c("GCM","rcp","proj.year"), sep = "_", remove = T)
Y0$proj.year <- gsub(".gcm","",Y0$proj.year)

Y0 <- addVars(Y0)

Y0 <- Y0 %>% dplyr::select(GCM, rcp,proj.year, vars)

##=======================
##Option 2 if folder with individual GCM outputs

#for(GCM in GCMs){

#Y0 <- fread(paste0("./inputs/",grid,"_",GCM,".csv", sep = ""), select = c("Year", "ID1","ID2", Columns), data.table = F)   ##
# Y0 <- separate(Y0, Year, into = c("GCM","rcp","proj.year"), sep = "_", remove = T)
# Y0$proj.year <- gsub(".gcm","",Y0$proj.year)
# 
# Y0 <- addVars(Y0)
# 
# Y0 <- Y0 %>% select(GCM, rcp,proj.year, vars)
##======================
require(doParallel)
set.seed(123321)
coreNum <- as.numeric(detectCores()-2)
cl <- makeCluster(coreNum)
registerDoParallel(cl, cores = coreNum)

out <- foreach(rcp = rcps, .combine = rbind) %:%
  foreach(proj.year = proj.years, .combine = rbind, .packages = c("data.table","ranger", "dplyr")) %dopar% {
    sub <- Y0[which(Y0$GCM == GCM & Y0$rcp == rcp & Y0$proj.year == proj.year),]
    subPred <- predict(BGCmodel, sub)
    
    fwrite(list(subPred$predictions), paste("outputs/BGC.pred", grid, GCM, rcp, proj.year, ".csv", sep = ""))
    # temp <- aggregate(ID ~ BGC.pred, data = subPred, FUN = length) %>% 
    #   mutate(GCM = GCM, rcp = rcp, proj.year = proj.year)
    # temp
    print(GCM)
  }

stopCluster(cl)

