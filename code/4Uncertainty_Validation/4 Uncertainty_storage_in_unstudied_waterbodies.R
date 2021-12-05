###################################################################################
##Estimating storage trend uncertainty in unstudied water bodies###################
###################################################################################
library(mblm)
library(stringr)
library(ggplot2)
library(mosaic)
library(Kendall)
options("scipen"=100, "digits"=16)

cat("For lakes with an area of 10-100 km2, using LakeID; For lakes smaller than 10 km2, using LakeIDold")
##input: Folder that contains preprocessed yearly area time series
inputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/Uncertainty/Eachlake/GT10km2/"
##input: preprocessed area variation (max.area-min.area) for each water body
Lake_arearangecsv = "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/Uncertainty/Glakes_GT10km2_JRCMin_Maxareas_inYears2019and2020.csv"
##input: preprocessed level variation (max.level-min.level) for each water body
ICESat2Levelragecsv="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/LakeLevel/ICESat2_Cooley/ForUncertainty/ICESat2_lake_variability_v1_Nov2020_wLakeID_Levelrange_UnStudiedLakes_GT10km2.csv"
##input: lake information csv 
LakeLookUpTablecsv="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/LakeLevel/ICESat2_Cooley/ForUncertainty/GlakesROI_GT10_points.csv"
##output: LWS rate for each water body
output="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/Uncertainty/0LWSRate_Uncertainty_GT10km2_LakeID.txt"

##Check invalid level variation
levelrange.invalidthres=200
#minimum points for Man Kendall trends 
minimumpoints.MKtrend=6

StartYear=1992
EndYear=2020

suppressPackageStartupMessages(source("C:/0UCBShare/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Lib_SenSlope.R"))
source("C:/0UCBShare/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Lib_SenSlope.R")

#extract number from a string
numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
} 

fileList <- list.files(inputDir,pattern="\\.csv$")

Levelrange.df <- data.frame()
Arearange <- read.csv(file = Lake_arearangecsv, header = TRUE, sep = ",")
LakeLookUpTable= read.csv(file = LakeLookUpTablecsv, header = TRUE, sep = ",")
ICESat2Levelrage <- read.csv(file = ICESat2Levelragecsv, header = TRUE, sep = ",")
ICESat2Levelrage$IsReservoir=ICESat2Levelrage$TypeID==3

ICESat2Levelrage.group <- 
    ICESat2Levelrage[ICESat2Levelrage$HYBAS_ID>0,] %>%
    group_by(SizeGroup,HYBAS_ID,IsReservoir) %>%
    summarise(Levelrangemedium = median(Levelrange))

sink(output)
cat("LakeIDandName LWSRate Pvalue Levelrange LevelrangeLT SizeGroup HYBAS_ID IsReservoir")
cat("\n")
for (inputfile in fileList) {
    #cat(inputfile)
    inputfilefull = paste(inputDir, inputfile, sep = "")
    poss=which(strsplit(inputfile, "")[[1]]=="_")
    subname=substr(inputfile, 0, poss[1]-1)
    cat(subname)
    cat(" ")
    LakeID <- numextract(substr(subname,3, nchar(subname)))
    
    subname = substr(inputfile, 0, nchar(inputfile) - 4)
 
    mydata <- read.csv(file = inputfilefull, header = TRUE, sep = ",")
    # Get all available Level and area data
    Area <- mydata$Area
    Year <- mydata$Year
    if(length(Area[!is.na(Area)])<2){
        next
    }
    Arearangetemp <- Arearange[Arearange$LakeID==LakeID,]
    if(length(Arearangetemp)==0){
        next
    }
    Areamax.seasonal <- Arearangetemp$MaxArea
    Areamin.seasonal <- Arearangetemp$MinArea
    ##Check invalid data in area
    if((Areamax.seasonal-Areamin.seasonal)<1e-3){
        next
    }
    
    ICESat2Levelragetemp <- ICESat2Levelrage[ICESat2Levelrage$LakeID==LakeID,]
    if(dim(ICESat2Levelragetemp)[1]==0){
        LakeLookUpTabletemp=LakeLookUpTable[LakeLookUpTable$LakeID==LakeID,]
        IsReservoirtemp = LakeLookUpTabletemp$TypeID==3
        ICESat2Levelrage.grouptemp=ICESat2Levelrage.group[ICESat2Levelrage.group$SizeGroup==LakeLookUpTabletemp$SizeGroup & ICESat2Levelrage.group$HYBAS_ID == 
                                                              LakeLookUpTabletemp$HYBAS_ID & ICESat2Levelrage.group$IsReservoir==IsReservoirtemp,]
        if(dim(ICESat2Levelrage.grouptemp)[1]==0){
            next
        }
        Levelrange=ICESat2Levelrage.grouptemp$Levelrangemedium
        ##Obtain info of lake size, basin, and lake type (reservoir or natural)
        SizeGroup=ICESat2Levelrage.grouptemp$SizeGroup
        HYBAS_ID=ICESat2Levelrage.grouptemp$HYBAS_ID
        IsReservoir=ICESat2Levelrage.grouptemp$IsReservoir
    }else{
        Levelrange <- ICESat2Levelragetemp$Levelrange
        SizeGroup=ICESat2Levelragetemp$SizeGroup
        HYBAS_ID=ICESat2Levelragetemp$HYBAS_ID
        IsReservoir=ICESat2Levelragetemp$IsReservoir
    }
    
    ##Check invalid data in area
    if((Areamax.seasonal-Areamin.seasonal)<1e-3 | Levelrange<1e-3){
        next
    }
    
    temp.df <- data.frame(LakeID=LakeID,Levelrange=Levelrange)
    Levelrange.df <- rbind(Levelrange.df, temp.df)
    
    kEval1 <- (Areamax.seasonal-Areamin.seasonal)/Levelrange 
    Intercept <- 0 
    ##Check if negative correlation between area and level
    isNslope=kEval1 <0

    f1 <- function(x) {
        Intercept + kEval1 * x 
    }
    f=f1
    
    YearlyArea.max=max(Area,na.rm = TRUE)
    YearlyArea.min=min(Area,na.rm = TRUE)
    Levelrange.LT = (YearlyArea.max-YearlyArea.min)/kEval1
    if(Levelrange.LT>levelrange.invalidthres){
        next
    }
    ##set possible lake level range (meter)
    ##This is relative level, change the limits if necessary
    uplimit=8848*2
    lowlimit=-1000
   
    r_square<--99
    adjr_square=r_square
    isNslope=FALSE
    if(isNslope){
       r_square=0
       adjr_square=0
    }

    subname=paste(subname,round(r_square,3),sep=" (R^2:")
    subname=paste(subname,")",sep="")

    #mse <- mean(residuals(fit)^2)
    rmse <- -99 #sqrt(mse)
    #rmse
    levelerrorOri = -99 #mydata$LevelErr
    #levelerror=levelerrorOri[!is.na(levelerrorOri)]
    funSq <- function(x) {
        x * x
    }
    #LerrorSq = funSq(levelerror)
    #LerrorMean = sqrt(sum(LerrorSq))
    LerrorMean = -99 #LerrorMean/length(LerrorSq)

    ##Get Inverse of the function f
    inverse = function(f, lower = -100, upper = 100) {
        function(y) uniroot((function(x) f(x) - y), lower = lower, 
            upper = upper)[1]
    }
    f_inverse = inverse(f, lowlimit,uplimit)

    lastlevel = -101
    lastV = 0
    isFirstvalid = 0
    
    vec_volume<-c()
    vec_year<-c()
    for (i in 1:length(Area)) {
        if (!is.na(Area[i])) {
            if (i > 1 & isFirstvalid > 0) {
                if (lastlevel < f_inverse(Area[i])$root) 
                    tempV = integrate(f, lastlevel, f_inverse(Area[i])$root)$val else tempV = -1 * (integrate(f, f_inverse(Area[i])$root, 
                                                                                                              lastlevel)$val)
                    lastV = tempV/1000 + lastV
                    vec_volume<-c(vec_volume,lastV)
                    vec_year<-c(vec_year,StartYear+i-1)
            } else {
                isFirstvalid = 1
                vec_volume<-c(vec_volume,lastV)
                vec_year<-c(vec_year,StartYear)
            }
            lastlevel = f_inverse(Area[i])$root
        } 
    }
    ## create a time series data type for the trend calculation
    volume.initialdata <- array(data = NA, dim = EndYear-StartYear+1, dimname = NULL)
    volume.ts <- ts(volume.initialdata,frequency = 1, start = c(StartYear))
    
    volume.ts[vec_year-StartYear+1]=vec_volume

    if(length(vec_volume)<minimumpoints.MKtrend)
    {
        fit <- lm(vec_volume~vec_year)
        xEval <- coef(summary(fit))["vec_year","Estimate"]
        xPval <- coef(summary(fit))["vec_year","Pr(>|t|)"]
    }else{
        xEval <- senSlope(vec_volume~vec_year)$coefficients[2]
        xPval <- MannKendall(volume.ts)$sl[1]
    }

    if(xPval<1e-4){
        xPval=0
    }
    cat(format(xEval, 
               scientific=F))
    cat(" ")
    cat(format(xPval, 
               scientific=F))
    cat(" ")
    cat(Levelrange)
    cat(" ")
    cat(Levelrange.LT)
    cat(" ")
    cat(SizeGroup)
    cat(" ")
    cat(format(HYBAS_ID, 
               scientific=F))
    cat(" ")
    cat(as.integer(IsReservoir))
    cat(" ")
    cat("\n")
}
sink()


