###################################################################################
##Estimating lake water storage using water area time series and ICESat-1 levels###
###################################################################################
library(mblm)
library(stringr)
library(ggplot2)
options("scipen"=100, "digits"=16)

##input: folder that contains preprocessed area and level time series
inputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/AandL/0AllLakes4Storage/LakeAreaandLevel_used/Natural/woUsableLevels/0excludedinCooley/0TryICESat1/"
##output: estimated lake water storage time series
outputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/GlakesAllTypes/WaterVolume/0Natural/"
fileList <- list.files(inputDir,pattern="\\_monthly.csv$")

cat("A-L match period: Jan 03 to Dec 09")
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
    outputfilefull = paste(outputfilefull, "Volume_ICESatLevelrange.csv", sep = "")
    
    mydata <- read.csv(file = inputfilefull, header = TRUE, sep = ",")
    # Get all available Level and area data
    #Timebg <- mydata$Time[133:216]
    Level <- mydata$LakeLevel
    Area <- mydata$LakeArea
    ## Non-Cloud cover
    NonCCover<-mydata$NonCCover
    ##Original area time series before cleaned for errors
    AreaOri<- mydata$LakeAreaOri
    Area.ICESat <- Area[133:216]
    if(length(Area.ICESat[!is.na(Area.ICESat)])<Numareas_Thres.hypsometryCali){
        next
    }
    Area.ICESatmax <- max(Area.ICESat,na.rm = TRUE)
    Area.ICESatmin <- min(Area.ICESat,na.rm = TRUE)
    
    
    
    Level.ICESat <- Level[133:216]
    if(length(Level.ICESat[!is.na(Level.ICESat)])==0){
        next
    }
    Level.ICESatmax <- max(Level.ICESat,na.rm = TRUE)
    Level.ICESatmin <- min(Level.ICESat,na.rm = TRUE)
    ##Check if no valid data available
    if((Area.ICESatmax-Area.ICESatmin)<1e-3 | (Level.ICESatmax-Level.ICESatmin)<1e-3){
        next
    }
    Levelrange=    Level.ICESatmax-Level.ICESatmin
    
    temp.df <- data.frame(LakeID=LakeID,Levelrange=Levelrange)
    Levelrange.df <- rbind(Levelrange.df, temp.df)
    
    kEval1 <- (Area.ICESatmax-Area.ICESatmin)/Levelrange 
    Intercept <- 0 
    ##Check if negative correlation between area and level
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

data.Levelrangecheck <- data.frame(
    name=c(rep("Level range",dim(Levelrange.df)[1])),
    value=c(Levelrange.df$Levelrange)
)
###Check if any evident errors
# data.Levelrangecheck %>%
#     ggplot( aes(x=name, y=value, fill=name)) +
#     geom_boxplot() +
#     scale_y_continuous(expand = c(0,0),limits = c(0,5))+ #breaks=c(0,0.2,0.4,0.6,0.8,1)
#     theme(
#         legend.position="none",
#         plot.title = element_text(size=9),
#         axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
#     ) +
#     ggtitle("") +
#     ylab("Trend (%/yr)") +
#     xlab("")
#Check for natural lakes; 10 m check for reservoirs
Levelrange.df.filter <- Levelrange.df[Levelrange.df$Levelrange>5,]

