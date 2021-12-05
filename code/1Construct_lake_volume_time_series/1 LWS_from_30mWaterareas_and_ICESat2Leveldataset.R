######################################################################################################################
##Estimating lake water storage using water area time series and ICESat-2 levels (from Cooley et al. 2021, Nature)####
######################################################################################################################
library(mblm)
library(stringr)
library(ggplot2)
options("scipen"=100, "digits"=16)

##input: folder that contains preprocessed area time series
inputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/AandL/0AllLakes4Storage/LakeAreaandLevel_used/Reservoirs/U211120/ICESat2/"
##input: level variation obtained from Cooley et al. 2021
ICESat2Levelragecsv="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/LakeLevel/ICESat2_Cooley/ICESat2_lake_height_time_series_v1_Nov2020_Levelrange_StudyLakes.csv"
##output: estimated lake water storage time series
outputDir="C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/GlakesAllTypes/WaterVolume/0Reservoirs/U211120/"

fileList <- list.files(inputDir, '_monthly.csv$', full.name = F)


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
    poss=which(strsplit(inputfile, "")[[1]]=="_")
    subname=substr(inputfile, 0, poss[1]-1)
    LakeID <- numextract(substr(subname,3, nchar(subname)))
    
    subname = substr(inputfile, 0, nchar(inputfile) - 4)
    outputfilefull = paste(outputDir, subname, sep = "")
    outputfilefull = paste(outputfilefull, "Volume_ICESat2Levelrange.csv", sep = "")
 
    mydata <- read.csv(file = inputfilefull, header = TRUE, sep = ",")

    Level <- mydata$LakeLevel
    Area <- mydata$LakeArea
    ## Non-contaminated proportion of satellite images
    NonCCover<-mydata$NonCCover
    ##Original area time series before filtering evident errors
    AreaOri<- mydata$LakeAreaOri
    ##Get areas from Oct 18 to Jul 20
    Area.ICESat2 <- Area[323:343]
    if(length(Area.ICESat2[!is.na(Area.ICESat2)])<Numareas_Thres.hypsometryCali){
        ##a small fraction of lakes have a long frozen season (>10 months), thus considering a longer period for A-L match: Apr 17 to Jul20 [304:343]
        Area.ICESat2 <- Area[304:343]
    }
    if(length(Area.ICESat2[!is.na(Area.ICESat2)])<Numareas_Thres.hypsometryCali){
        next
    }
    Area.ICESat2max <- max(Area.ICESat2,na.rm = TRUE)
    Area.ICESat2min <- min(Area.ICESat2,na.rm = TRUE)
    
    ICESat2Levelrange <- read.csv(file = ICESat2Levelragecsv, header = TRUE, sep = ",")
    ICESat2Levelrange <- ICESat2Levelrange[ICESat2Levelrange$LakeID==LakeID,]
    ##Check if no valid data available
    if(dim(ICESat2Levelrange)[1]==0 | (Area.ICESat2max-Area.ICESat2min)<1e-3){
        next
    }

    temp.df <- data.frame(LakeID=LakeID,Levelrange=ICESat2Levelrange$Levelrange)
    Levelrange.df <- rbind(Levelrange.df, temp.df)
    
    kEval1 <- (Area.ICESat2max-Area.ICESat2min)/ICESat2Levelrange$Levelrange 
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

    rmse <- -99 #sqrt(mse)
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
    ##Only for the study period 199210-202009
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
###Check if any evident errors
data.Levelrangecheck <- data.frame(
    name=c(rep("Level range",dim(Levelrange.df)[1])),
    value=c(Levelrange.df$Levelrange)
)
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
#Check for natural lakes; 10 m if check for reservoirs
Levelrange.df.filter <- Levelrange.df[Levelrange.df$Levelrange>5,]

