###################################################################################
##Calculating lake water storage trend using Mann-Kendall method###################
###################################################################################
setwd("G:/My Drive/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Glakes/")

library(XLConnect)
library(mblm)
require(ggplot2)
library(Kendall)

options("scipen"=100, "digits"=8)

##input: folders that contains lake volume time series
inputDir="./data/Output_Volume_TimeSeries/Monthly/"
##input: Time steps to be studied
TimeTemplateCsv="./data/Preprocessed/Template/0TimeTemplateforLandsatMonthly_199210_202009.csv"
##output: LWS trends
outputtxt='./data/Output_Volume_TimeSeries/Output_lake_volume_trends_1992_2020.txt'

##Refer to https://rdrr.io/github/USGS-R/smwrStats/man/senSlope.html
suppressPackageStartupMessages(source("C:/0UCBShare/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Lib_SenSlope.R"))
source("C:/0UCBShare/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Lib_SenSlope.R")

###########Read the time steps##############
myTimedata <- read.csv(file=TimeTemplateCsv, header=TRUE, sep=",")
Time=myTimedata$Time   
## as continuous integer with step of 1
TimeNum=myTimedata$TimeNum 

fileList <- list.files(inputDir, '.csv$', full.name = F)
##Write to file
sink(outputtxt)
cat("LakeIDandName LWSRate Pvalue MonthlyCoverage Rsquare RMSE LevelError TrendR2 MaximumValue")
cat("\n")
for(inputfile in fileList)
{
  inputfilefull=paste(inputDir,inputfile,sep="")
  ## find the first "_", bef is the lake Name
  poss=regexpr("30m",inputfile)[1]
  subname=substr(inputfile, 0, poss[1]-1)
  cat(subname)
  cat(" ")
  ##read the first line
  con <- file(inputfilefull, "r")
  firstline=readLines(con,n=1)
  close(con)
  
  ## Get other information of LWS data
  posRMSE=regexpr("RMSE of relation:",firstline)[1]
  RMSE=substring(firstline,posRMSE+17)
  RMSE=substring(RMSE,0,regexpr(",",RMSE)[1]-1)
  posLevelError=regexpr("LevelError",firstline)[1]
  LevelError=substring(firstline,posLevelError+11)
  LevelError=substring(LevelError,0,regexpr(",",LevelError)[1]-1)
  posRsquare=regexpr("Multiple R-squared:",firstline)[1]
  Rsquare=substring(firstline,posRsquare+19)
  Rsquare=substring(Rsquare,0,regexpr(",",Rsquare)[1]-1)
  if(regexpr("HModel",inputfile)[1]>-1){
    RMSE=-99
    Rsquare=-99
    LevelError=-99
  }
  
  mydata <- read.csv(file=inputfilefull, header=TRUE, sep=",")
  
  ## create a time series data type for the trend calculation
  volume.initialdata <- array(data = NA, dim = length(TimeNum), dimname = NULL)
  volume.ts <- ts(volume.initialdata,frequency = 12, start = c(1992, 10))
  ##LWS presents water storage here    
  LWS=mydata$rws
  LWS=LWS[1:length(LWS)]
  LWSsel=LWS[!is.na(LWS)]
  Timesel=Time[!is.na(LWS)]
  TimeNumsel=TimeNum[!is.na(LWS)]
  
  ## skip NA values
  LWSinDeed=mydata$LakeLWS
  LWSinDeed=LWSinDeed[!is.na(LWSinDeed)]
  
  MonthlyCoverage=length(Timesel)/length(Time)
  
  ##Checking for potential errors (e.g., constant LWS values) 
  if((max(LWSinDeed)-min(LWSinDeed)<1e-3)&((max(LWSsel)-min(LWSsel)<1e-3))){
    cat(0)
    cat(" ")
    cat(1)
    cat(" ")
    cat(formatC(MonthlyCoverage, 
                digits = 3))
    cat(" ")
    cat(format(Rsquare, 
               scientific=F))
    cat(" ")
    cat(format(RMSE, 
               scientific=F))
    cat(" ")
    cat(format(LevelError, 
               scientific=F))
    cat(" ")
    cat(0)
    cat(" ")
    if(length(LWSinDeed)==0)
    {
      cat("-9999")
    }
    else
    {
      cat(max(LWSinDeed))
    }
    cat("\n")
    next
  }
  
  LWSdata<-data.frame(LWS=LWSsel,time=as.Date(Timesel, format = "%m/%d/%Y"),timenum=TimeNumsel)
  ##Deseasonality by substracting the long-term mean in each month
  LWSdatamon01=subset(LWSdata, format.Date(time, "%m")=="01")
  deLWSmon01=LWSdatamon01$LWS-mean(LWSdatamon01$LWS)
  LWSdatamon01$LWS <- deLWSmon01
  LWSdatamon02=subset(LWSdata, format.Date(time, "%m")=="02")
  deLWSmon02=LWSdatamon02$LWS-mean(LWSdatamon02$LWS)
  LWSdatamon02$LWS <- deLWSmon02
  LWSdatamon03=subset(LWSdata, format.Date(time, "%m")=="03")
  deLWSmon03=LWSdatamon03$LWS-mean(LWSdatamon03$LWS)
  LWSdatamon03$LWS <- deLWSmon03
  LWSdatamon04=subset(LWSdata, format.Date(time, "%m")=="04")
  deLWSmon04=LWSdatamon04$LWS-mean(LWSdatamon04$LWS)
  LWSdatamon04$LWS <- deLWSmon04
  LWSdatamon05=subset(LWSdata, format.Date(time, "%m")=="05")
  deLWSmon05=LWSdatamon05$LWS-mean(LWSdatamon05$LWS)
  LWSdatamon05$LWS <- deLWSmon05
  LWSdatamon06=subset(LWSdata, format.Date(time, "%m")=="06")
  deLWSmon06=LWSdatamon06$LWS-mean(LWSdatamon06$LWS)
  LWSdatamon06$LWS <- deLWSmon06
  LWSdatamon07=subset(LWSdata, format.Date(time, "%m")=="07")
  deLWSmon07=LWSdatamon07$LWS-mean(LWSdatamon07$LWS)
  LWSdatamon07$LWS <- deLWSmon07
  LWSdatamon08=subset(LWSdata, format.Date(time, "%m")=="08")
  deLWSmon08=LWSdatamon08$LWS-mean(LWSdatamon08$LWS)
  LWSdatamon08$LWS <- deLWSmon08
  LWSdatamon09=subset(LWSdata, format.Date(time, "%m")=="09")
  deLWSmon09=LWSdatamon09$LWS-mean(LWSdatamon09$LWS)
  LWSdatamon09$LWS <- deLWSmon09
  LWSdatamon10=subset(LWSdata, format.Date(time, "%m")=="10")
  deLWSmon10=LWSdatamon10$LWS-mean(LWSdatamon10$LWS)
  LWSdatamon10$LWS <- deLWSmon10
  LWSdatamon11=subset(LWSdata, format.Date(time, "%m")=="11")
  deLWSmon11=LWSdatamon11$LWS-mean(LWSdatamon11$LWS)
  LWSdatamon11$LWS <- deLWSmon11
  LWSdatamon12=subset(LWSdata, format.Date(time, "%m")=="12")
  deLWSmon12=LWSdatamon12$LWS-mean(LWSdatamon12$LWS)
  LWSdatamon12$LWS <- deLWSmon12
  LWSdatamonAll=Reduce(function(x, y) merge(x, y, all=TRUE), list(LWSdatamon01,LWSdatamon02,LWSdatamon03,LWSdatamon04,LWSdatamon05, LWSdatamon06, LWSdatamon07,LWSdatamon08,LWSdatamon09, LWSdatamon10, LWSdatamon11,LWSdatamon12))
  #LWSdatamonAll[order(LWSdatamonAll$timenum),]
  DeLWSsel=LWSdatamonAll$LWS
  Detimenumsel=LWSdatamonAll$timenum
  ##using lm for trend only for comparing
  fit <- lm(DeLWSsel~Detimenumsel)
  #using Mannkendall for trend
  volume.ts[Detimenumsel]=DeLWSsel
  ###mk.test in the Trend package requires complete time series; while Mankendall in the Kendall package can allow missing values
  #MannKendall(volume.ts)$sl[1]
  #sens.slope requires complete time series; while senSlope (developed by USGS Water resources) does not:https://rdrr.io/github/USGS-R/smwrStats/man/senSlope.html
  #senSlope(DeLWSsel~Detimenumsel)$coefficients[2]
  
  #plot(volume.ts)
  #plot(Detimenumsel,DeLWSsel)
  
  xEval <- senSlope(DeLWSsel~Detimenumsel)$coefficients[2]
  xPval <- MannKendall(volume.ts)$sl[1]

  if(xPval<1e-4){
    xPval=0
  }
  cat(format(xEval*12, 
             scientific=F))
  cat(" ")
  cat(format(xPval, 
             scientific=F))
  cat(" ")
  cat(formatC(MonthlyCoverage, 
              digits = 3))
  cat(" ")
  cat(format(Rsquare, 
             scientific=F))
  cat(" ")
  cat(format(RMSE, 
             scientific=F))
  cat(" ")
  cat(format(LevelError, 
             scientific=F))
  cat(" ")
  cat(format(summary(fit)$r.squared, 
             scientific=F))
  cat(" ")
  if(length(LWSinDeed)==0)
  {
    cat("-9999")
  }
  else
  {
    cat(max(LWSinDeed))
  }
  cat("\n")
  
}
sink()

