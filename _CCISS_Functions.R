##Functions called from various CCISS scripts

addVars <- function(dat){
  dat$PPT_MJ <- dat$PPT05 + dat$PPT06  # MaY/June precip
  dat$PPT_JAS <- dat$PPT07 + dat$PPT08 + dat$PPT09  # July/Aug/Sept precip
  dat$PPT.dormant <- dat$PPT_at + dat$PPT_wt  # for calculating spring deficit
  dat$CMD.def <- 500 - (dat$PPT.dormant)  # start of growing season deficit original value was 400 but 500 seems better
  dat$CMD.def[dat$CMD.def < 0] <- 0  #negative values set to zero = no deficit
  dat$CMDMax <- dat$CMD07
  dat$CMD.total <- dat$CMD.def + dat$CMD
  return(dat)
}