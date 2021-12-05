###################################################################################
##Estimating lake water storage using water area time series and Hypsometry data ##
##from Li et al. 2020, Remote Sensing of Environment###############################
###################################################################################
library(mblm)
library(stringr)
options("scipen"=100, "digits"=16)

inputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/AandL/0AllLakes4Storage/LakeAreaandLevel_used/Reservoirs/U211120/Gao/"
outputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/GlakesAllTypes/WaterVolume/0Reservoirs/U211120/"
fileList <- list.files(inputDir,pattern="\\_monthly.csv$")
ReservoirHypsometrybyGao="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/LakeLevel/Bathymetry_Gao/0_Global_reservoir_bathymetry_dataset_onStudyLakes.csv"

##extract number from a string
numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
} 

for (inputfile in fileList) {
    cat(inputfile)
    inputfilefull = paste(inputDir, inputfile, sep = "")
    poss=which(strsplit(inputfile, "")[[1]]=="_")
    subname=substr(inputfile, 0, poss[1]-1)
    LakeID <- numextract(substr(subname,3, nchar(subname)))
    
    subname = substr(inputfile, 0, nchar(inputfile) - 4)
    outputfilefull = paste(outputDir, subname, sep = "")
    outputfilefull = paste(outputfilefull, "Volume_LinearGao.csv", sep = "")
 
    mydata <- read.csv(file = inputfilefull, header = TRUE, sep = ",")
    ## Get all available Level and area data
    Level <- mydata$LakeLevel
    Area <- mydata$LakeArea
    ## Non-Cloud cover
    NonCCover<-mydata$NonCCover
    ##Area time series before cleaned for errors
    AreaOri<- mydata$LakeAreaOri
    
    myhypsometry <- read.csv(file = ReservoirHypsometrybyGao, header = TRUE, sep = ",")
    myhypsometry <- myhypsometry[myhypsometry$LakeID==LakeID,]
    if(dim(myhypsometry)[1]==0){
        next
    }
    ##set possible lake level range (meter)
    minL=-1000
    maxL=8848
    
    slope.Gao <- myhypsometry$a
    intercept.Gao <- myhypsometry$b
    
    kEval1 <- 1.0/slope.Gao 
    Intercept <- -1.0*intercept.Gao/slope.Gao
    ##Check if negative correlation between area and level
    isNslope=kEval1 <0

    f1 <- function(x) {
        Intercept + kEval1 * x 
    }
    f=f1
    ##set possible lake level range (meter)
    uplimit=8848
    lowlimit=-1000
   
    r_square<-myhypsometry$r2
    adjr_square=r_square
    isNslope=FALSE
    if(isNslope){
       r_square=0
       adjr_square=0
    }

    subname=paste(subname,round(r_square,3),sep=" (R^2:")
    subname=paste(subname,")",sep="")
    #plot(LakeLevelsel,LakeAreasel,main=subname,xlab="Lake Level (m)",ylab="Lake Area (km^2)")
    #lines(Seq_Level,f(Seq_Level))
    summary(fit)$r.squared
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

    ##Output the lake volume time series
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
    #############################################################################################################
    ################19921001-20200901
    date <- mydata$Time[10:345]
    levelerrorOri=levelerrorOri[10:345]
    Level=Level[10:345]
    Area=Area[10:345]
    AreaOri=AreaOri[10:345]
    NonCCover=NonCCover[10:345]
    #NonCCover=NonCCover[10:322]
    # 
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

