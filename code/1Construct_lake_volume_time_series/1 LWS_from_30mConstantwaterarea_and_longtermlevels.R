########################################################################################
##Estimating lake water storage using long-term level time series and a constant area###
########################################################################################
library(XLConnect)
library(mblm)
options("scipen"=100, "digits"=16)

##input: folder that contains preprocessed area and level time series
inputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/AandL/0AllLakes4Storage/LakeAreaandLevel_used/Reservoirs/wUsableLevels/LevelDominant/LT060/"
##output: estimated lake water storage time series
outputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/GlakesAllTypes/WaterVolume/0Reservoirs/"

fileList <- list.files(inputDir,pattern="\\_monthly.csv$")#_monthly.csv
for (inputfile in fileList) {
    cat(inputfile)
    inputfilefull = paste(inputDir, inputfile, sep = "")
    subname = substr(inputfile, 0, nchar(inputfile) - 4)
    outputfilefull = paste(outputDir, subname, sep = "")
    outputfilefull = paste(outputfilefull, "Volume_ConstantArea.csv", sep = "")

    mydata <- read.csv(file = inputfilefull, header = TRUE, sep = ",")
    # Get all available Level and area data
    Level <- mydata$LakeLevel
    Area <- mydata$LakeArea
    ##Only for the study period 199210-202009
    Areasub<-Area [10:345]
    MeanA=mean(Areasub[!is.na(Areasub)])
    ## Non-contaminated proportion of satellite images
    NonCCover<-mydata$NonCCover
    ##Original area time series before filtering evident errors
    AreaOri<- mydata$LakeAreaOri
    ##Areas used for calibrating hypsometry: areas from good images only or areas from good and contaminated (bad) images
    selareafield=mydata$AreaInHypsometry[1]
    if (selareafield[1]=="goodonly") {
        LakeAreatemp=Area
        LakeAreatemp[which(NonCCover < 0.95, arr.ind = TRUE)]<- NA
        LakeLeveltemp=Level
        LakeLeveltemp[which(NonCCover < 0.95, arr.ind = TRUE)]<- NA
        LakeAreasel=LakeAreatemp[!is.na(LakeAreatemp)&!is.na(LakeLeveltemp)]
        LakeLevelsel=LakeLeveltemp[!is.na(LakeAreatemp)&!is.na(LakeLeveltemp)]
    } else if (selareafield[1]=="good_bad"){
        LakeAreatemp=Area
        LakeLeveltemp=Level
        LakeAreasel=Area[!is.na(LakeAreatemp)&!is.na(LakeLeveltemp)]
        LakeLevelsel=Level[!is.na(LakeAreatemp)&!is.na(LakeLeveltemp)]
    }
    
    minL=min(Level[!is.na(Level)])
    maxL=max(Level[!is.na(Level)])

    ##get areas for filling the gaps of levels and subsequently volumes
    Areaforvolume=Area[is.na(Level)&!is.na(Area)]
    ##set possible lake level range (meter)
    uplimit=8848
    lowlimit=-1000
    ##for checking inflection points of the hypsometric curve
    Lrange<-seq(from=lowlimit, to=uplimit, by=0.1)
    r_square=0
    adjr_square=0
    rmse =0
    if(length(LakeAreasel)>2){
       fit1<- lm(LakeAreasel~LakeLevelsel) 
       kEval1 <- coef(summary(fit1))["LakeLevelsel", "Estimate"]
       Intercept <- coef(summary(fit1))["(Intercept)", "Estimate"]
       ##Check if negative correlation between area and level
       isNslope=kEval1 <0

       f1 <- function(x) {
           Intercept + kEval1 * x 
       }
       f=f1
       fit=fit1
       Seq_Level <- seq(from=minL, to=maxL, by=0.1)

       r_square=summary(fit)$r.squared
       adjr_square=summary(fit)$adj.r.squared
       if(isNslope){
          r_square=0
          adjr_square=0
       }
       mse <- mean(residuals(fit)^2)
       rmse <- sqrt(mse)
    }
    
    subname=paste(subname,round(r_square,3),sep=" (R^2:")
    subname=paste(subname,")",sep="")

    levelerrorOri = mydata$LevelErr
    levelerror=levelerrorOri[!is.na(levelerrorOri)]
    funSq <- function(x) {
        x * x
    }
    LerrorSq = funSq(levelerror)
    LerrorMean = sqrt(sum(LerrorSq))
    LerrorMean = LerrorMean/length(LerrorSq)
    LerrorMean
    
    ##Output the lake volume time series
    sink(outputfilefull)
    cat("Time,Mean_Levels,Mean_Level_Errors,LakeArea,NonCCover,Cleaned_Area,wsc,rws,,,")
    cat("Area=")
    cat(MeanA)
    cat(",,,")
    cat("RMSE of relation: ", formatC(rmse, digits = 6), ",,,\tLevelError: ", 
        formatC(LerrorMean, digits = 6), ",,,", sep = " ")
    cat("Multiple R-squared: ", formatC(r_square, 
        digits = 3), ",,,\tAdjusted R-squared: ", formatC(adjr_square, 
        digits = 3), "\n", sep = " ")

    ##Only for the study period 199210-202009
    date <- mydata$Time[10:345]
    levelerrorOri=levelerrorOri[10:345]
    Level=Level[10:345]
    Area=Area[10:345]
    AreaOri=AreaOri[10:345]
    NonCCover=NonCCover[10:345]
    # 
    lastlevel = -101
    lastV = 0
    isFirstvalid = 0

    for (i in 1:length(date)) {
        cat(date[i])
        cat(",")
        if (!is.na(Level[i])) {
            cat(Level[i])
            cat(",")
            cat(levelerrorOri[i])
            cat(",")
            if(!is.na(AreaOri[i])){
                cat(AreaOri[i])
                cat(",")
                cat(NonCCover[i])
                if(!is.na(Area[i])){
                    cat(",")
                    cat(Area[i])
                }
                else {
                  cat(",")
                }
            }
            else {
              cat(",,")
            }
            if (i > 1 & isFirstvalid > 0) {
                cat(",")
                tempV=(Level[i]-lastlevel)*MeanA
                cat(tempV/1000)
                cat(",")
                cat(tempV/1000 + lastV)
                lastV = tempV/1000 + lastV
            } else {
                cat(",0,0")
                isFirstvalid = 1
            }
            lastlevel = Level[i]
        } else if (!is.na(Area[i])) {
            cat(",,")
            cat(AreaOri[i])
            cat(",")
            cat(NonCCover[i])
            cat(",")
            cat(Area[i])
            cat(",,")
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

