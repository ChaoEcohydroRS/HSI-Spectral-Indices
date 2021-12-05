########################################################################################################################
##Calculating aridity trend in major hydrological basins with largest lake volume changes###############################
##aridity is quantified using aridity index, a ratio of annual precipitation to annual potential evapotranspiration#####
##trend is calculated using the Mann-Kendall method ####################################################################
########################################################################################################################
library(ncdf4)
library(Thermimage)
library(raster)
library(doParallel)
library(foreach)
library(lubridate)
library(stringr)
library(Kendall)
options("scipen"=100, "digits"=8)

##input: root folder for data of hydroclimate and human variables for each hydrological basin
ClimateDataRootFolder <- 'C:/0D/Research/3GLakeWatch/AnnualClimate/ClimateOnHydroBasin/'
##output: basinwide trends in aridity
OutAItrend.txt <- 'C:/0D/Research/3GLakeWatch/AnnualClimate/ClimateOnHydroBasin/AridityTrend_ensemblePandPET_40basins.txt'

cat("Creat subfolders in the ClimateDataRootFolder and input datasets of yearly Prep and PET as csv files in each subfolder\n")
Prep.DS <- c("PtrMSWEP","PtrERA","PtrMERRA","PtrCRU") 
PET.DS <- c("PEtrCRU","PEtrGLEAM")

##Period for calculating aridity trends
BeginYear = 1992
EndYear = 2020


AnnualPrepcsvfolder <- paste(ClimateDataRootFolder,Prep.DS[1],sep="/")
fileList <- list.files(AnnualPrepcsvfolder, '.csv$', full.name = F)
Colnames_drop <- c("Year")

sink(OutAItrend.txt)
cat("LakeID AIRate Pvalue")
cat("\n")

for(inputfile in fileList)
{
  lakeid <- as.numeric(gsub(".*?([0-9]+).*", "\\1", inputfile))
  Pensemble_df <- data.frame(Year=seq(BeginYear, EndYear, by=1))
  PETensemble_df <- data.frame(Year=seq(BeginYear, EndYear, by=1))
  for(Prepdata in Prep.DS)
  {
    inputfilefull.prep=paste0(ClimateDataRootFolder,Prepdata,"/",inputfile)
    mydata.prep <- read.csv(file=inputfilefull.prep, header=TRUE, sep=",")
    Pensemble_df <- cbind(Pensemble_df,mydata.prep$Value)
  }
  for(PETdata in PET.DS)
  {
    inputfilefull.pet=paste0(ClimateDataRootFolder,PETdata,"/",inputfile)
    mydata.pet <- read.csv(file=inputfilefull.pet, header=TRUE, sep=",")
    PETensemble_df <- cbind(PETensemble_df,mydata.pet$Value)
  }
  Pensemble_df = Pensemble_df[,!(names(Pensemble_df) %in% Colnames_drop)]
  PETensemble_df = PETensemble_df[,!(names(PETensemble_df) %in% Colnames_drop)]
  Pensemble_df.mean = rowSums(Pensemble_df)/dim(Pensemble_df)[2]
  PETensemble_df.mean = rowSums(PETensemble_df)/dim(PETensemble_df)[2]
  
  AI_df = Pensemble_df.mean/PETensemble_df.mean
  
  AI.initialdata <- array(data = NA, dim = length(AI_df[1:23]), dimname = NULL)
  AI.ts <- ts(AI.initialdata,frequency = 1, start = c(BeginYear))
  
  AI.ts[mydata.prep$Year[1:23]-BeginYear+1]=AI_df[1:23]
  MannKendall(AI.ts)
  xEval <- MannKendall(AI.ts)$tau[1]
  xPval <- MannKendall(AI.ts)$sl[1]
  
  if(xPval<1e-4){
    xPval=0
  }
  cat(lakeid)
  cat(" ")
  cat(format(xEval, 
             scientific=F))
  cat(" ")
  cat(format(xPval, 
             scientific=F))
  cat(" ")
  cat("\n")
}
sink()

