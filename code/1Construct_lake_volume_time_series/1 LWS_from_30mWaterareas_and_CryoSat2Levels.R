###################################################################################
##Estimating lake water storage using water area time series and CryoSat-2 levels##
###################################################################################
library(mblm)
library(stringr)
library(ggplot2)
options("scipen"=100, "digits"=16)

##input: folder that contains preprocessed area time series
inputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/AandL/0AllLakes4Storage/LakeAreaandLevel_used/Reservoirs/U211120/CryoSat2/"
##input: preprocessed water level variation from CryoSat-2
CryoSat2Levelragecsv="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/LakeLevel/CryoSat2_Jeff/CryoSat2_lake_height_time_series_Levelrange_StudyLakes.csv"
##output: estimated lake water storage time series
outputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/GlakesAllTypes/WaterVolume/0Reservoirs/U211120/"

fileList <- list.files(inputDir, '_monthly.csv$', full.name = F)

cat("using Apr10-Sep20 for A-L match as CryoSat2 available from apr 2010")
##Minimum number of area observations required for hypsometry calibration
Numareas_Thres.hypsometryCali=3

#extract number from a string
numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
} 
Levelrange.df <- data.frame()

for (inputfile in fileList) {
    cat(inputfile)
    inputfilefull = paste(inputDir, inputfile, sep = "")
    ## find the first "_"
    poss=which(strsplit(inputfile, "")[[1]]=="_")
    subname=substr(inputfile, 0, poss[1]-1)
    LakeID <- numextract(substr(subname,3, nchar(subname)))
    
    subname = substr(inputfile, 0, nchar(inputfile) - 4)
    outputfilefull = paste(outputDir, subname, sep = "")
    outputfilefull = paste(outputfilefull, "Volume_CryoSat2Levelrange.csv", sep = "")
    
    mydata <- read.csv(file = inputfilefull, header = TRUE, sep = ",")
    ## Get all available Level and area data
    Level <- mydata$LakeLevel
    Area <- mydata$LakeArea
    ## Non-Cloud cover 
    NonCCover<-mydata$NonCCover
    ##Original area time series before cleaned for errors
    AreaOri<- mydata$LakeAreaOri

    Area.CryoSat2 <- Area[220:345]
    if(length(Area.CryoSat2[!is.na(Area.CryoSat2)])<Numareas_Thres.hypsometryCali){
        next
    }
    Area.CryoSat2max <- max(Area.CryoSat2,na.rm = TRUE)
    Area.CryoSat2min <- min(Area.CryoSat2,na.rm = TRUE)
    
    CryoSat2Levelrange <- read.csv(file = CryoSat2Levelragecsv, header = TRUE, sep = ",")
    CryoSat2Levelrange <- CryoSat2Levelrange[CryoSat2Levelrange$LakeID==LakeID,]

    if(dim(CryoSat2Levelrange)[1]==0 | (Area.CryoSat2max-Area.CryoSat2min)<1e-3 ){
        next
    }
    
    temp.df <- data.frame(LakeID=LakeID,Levelrange=CryoSat2Levelrange$Levelrange)
    Levelrange.df <- rbind(Levelrange.df, temp.df)
    #fit1<- lm(LakeAreasel~LakeLevelsel) 
    minL=-1000#min(Level[!is.na(Level)])
    maxL=8848#max(Level[!is.na(Level)])
    #get area without level
    
    
    kEval1 <- (Area.CryoSat2max-Area.CryoSat2min)/CryoSat2Levelrange$Levelrange #coef(summary(fit1))["LakeLevelsel", "Estimate"]
    Intercept <- 0 #coef(summary(fit1))["(Intercept)", "Estimate"]
    
    isNslope=kEval1 <0
    
    f1 <- function(x) {
        Intercept + kEval1 * x 
    }
    f=f1
    ##set possible lake level range (meter)
    uplimit=8848
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
    LerrorMean
    
    ##Get Inverse of the function f
    inverse = function(f, lower = -100, upper = 100) {
        function(y) uniroot((function(x) f(x) - y), lower = lower, 
                            upper = upper)[1]
    }
    f_inverse = inverse(f, lowlimit,uplimit)
    
    
    ##Output the lake volume time series
    sink(outputfilefull)
    cat("Time,Mean_Levels,Mean_Level_Errors,LakeArea,NonCCover,Cleaned_Area,wsc,rws,,,")
    cat("Area=")
    cat(Intercept)
    if(kEval1>0){ cat("+") }
    cat(kEval1)
    cat("*")
    cat("x")
    
    cat(",,,")
    cat("RMSE of relation: ", formatC(rmse, digits = 6), ",,,\tLevelError: ", 
        formatC(LerrorMean, digits = 6), ",,,", sep = " ")
    cat("Multiple R-squared: ", formatC(r_square, 
                                        digits = 3), ",,,\tAdjusted R-squared: ", formatC(r_square, 
                                                                                          digits = 3), "\n", sep = " ")
    ##Only for the study period 19921001-20200901
    date <- mydata$Time[10:345]
    levelerrorOri=levelerrorOri[10:345]
    Level=Level[10:345]
    Area=Area[10:345]
    AreaOri=AreaOri[10:345]
    NonCCover=NonCCover[10:345]

    lastlevel = -101
    lastV = 0
    isFirstvalid = 0
    for (i in 1:length(date)) {
        cat(date[i])
        cat(",")
        if (!is.na(Area[i])) {
            cat(f_inverse(Area[i])$root)
            cat(",,")
            cat(AreaOri[i])
            cat(",")
            cat(NonCCover[i])
            cat(",")
            cat(Area[i])
            if (i > 1 & isFirstvalid > 0) {
                cat(",")
                if (lastlevel < f_inverse(Area[i])$root) 
                    tempV = integrate(f, lastlevel, f_inverse(Area[i])$root)$val else tempV = -1 * (integrate(f, f_inverse(Area[i])$root, 
                                                                                                              lastlevel)$val)
                    cat(tempV/1000)
                    cat(",")
                    cat(tempV/1000 + lastV)
                    lastV = tempV/1000 + lastV
            } else {
                cat(",0,0")
                isFirstvalid = 1
            }
            lastlevel = f_inverse(Area[i])$root
        } else if (!is.na(AreaOri[i])) {
            cat(",,")
            cat(AreaOri[i])
            cat(",")
            cat(NonCCover[i])
            cat(",")
            cat(",")
            cat(",")
        }
        cat("\n")
    }
    sink()
}


