########################################################################################
##Estimating lake water storage using time-varying area and level observations##########
########################################################################################
setwd("C:/Workstation/PreviousProject/GlobalRivers/Glakes/")

library(XLConnect)
library(mblm)
library(gridExtra)
library(ggplot2)
options("scipen"=100, "digits"=16)

##input: folder that contains preprocessed area and level time series
inputDir="./data/Preprocessed/A_L_TimeSeries/"
##output: estimated lake water storage time series
outputDir="./data/Output_Volume_TimeSeries/Monthly/"

##Good fitting in hypsometry
R2inHypsometry_thres=0.6

##Scale km2 to m2
Scalekm2_m2=1e6

##Mean error in mapped lake areas: 2.2%, obtained from Yao et al. 2019, Remote Sensing of Environment
err.mappedarea = 0.022

fileList <- list.files(inputDir, 'monthly.csv$', full.name = FALSE)
for (inputfile in fileList) {
    cat(inputfile)
    inputfilefull = paste(inputDir, inputfile, sep = "")
    subname = substr(inputfile, 0, nchar(inputfile) - 4)
    outputfilefull = paste(outputDir, subname, sep = "")
    outputfilefull = paste(outputfilefull, "Volume_Poly.csv", sep = "")

    mydata <- read.csv(file = inputfilefull, header = TRUE, sep = ",")
    # Get all available Level and area data
    Level <- mydata$LakeLevel
    Level.err <- mydata$LevelErr
    p1 <- ggplot(mydata, aes(x=TimeNum/12+1992.833, y=Level)) + 
        ggtitle("Water levels from altimetry-derived water level databases") +
        xlab("Year") +
        ylab("Water level (m)") +
        geom_line() +
        geom_point()+
        geom_errorbar(aes(ymin=Level-Level.err, ymax=Level+Level.err), width=.2,position=position_dodge(0.05))
    
    p2 <- ggplot(mydata, aes(x=TimeNum/12+1992.833, y=LakeArea*Scalekm2_m2)) + 
        ggtitle("Water areas mapped from 30-m Landsat images") +
        xlab("Year") +
        ylab("Water area (m2)") +
        geom_line() +
        geom_point()+
        geom_errorbar(aes(ymin=LakeArea-LakeArea*err.mappedarea, ymax=LakeArea+LakeArea*err.mappedarea), width=.2,position=position_dodge(0.05))+
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

    
    
    Area <- mydata$LakeArea*Scalekm2_m2
    ## Non-contaminated proportion of satellite images
    NonCCover<-mydata$NonCCover
    ##Original area time series before filtering evident errors
    AreaOri<- mydata$LakeAreaOri*Scalekm2_m2
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
    ##linear fitting
    fit1<- lm(LakeAreasel~LakeLevelsel) 
    minL=min(Level[!is.na(Level)])
    maxL=max(Level[!is.na(Level)])
    ##get areas for filling the gaps of levels and subsequently volumes
    Areaforvolume=Area[is.na(Level)&!is.na(Area)]
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
    ##set possible lake level range (meter)
    uplimit=8848
    lowlimit=-1000
    ##for checking inflection points of the hypsometric curve
    Lrange<-seq(from=lowlimit, to=uplimit, by=0.1)
   
    kEval2=0
    kEval3=0
    ##Try Quadratic equation
    ##To avoid singularity issue, minus the minimum value
    Level1=LakeLevelsel-min(LakeLevelsel)
    Level2=LakeLevelsel^2-min(LakeLevelsel^2)
    Level3=LakeLevelsel^3-min(LakeLevelsel^3)
    fit2shift<-lm(LakeAreasel~Level1+Level2) 
    ##Using AIC to avoid overfitting
    if(length(coef(summary(fit2shift)))>8 & AIC(fit2shift)<AIC(fit) & !(isNslope)){
        kEval2_1temp <- coef(summary(fit2shift))["Level1", "Estimate"]    
        kEval2_2temp <- coef(summary(fit2shift))["Level2", "Estimate"] 
        Intercept2_temp <- coef(summary(fit2shift))["(Intercept)", "Estimate"]
        Intercept2_temp=Intercept2_temp-kEval2_1temp*min(LakeLevelsel)-kEval2_2temp*min(LakeLevelsel^2)
        if(!is.na(kEval2_2temp)&(summary(fit2shift)$r.squared>R2inHypsometry_thres))
        {
            f2 <- function(x) {
            Intercept2_temp + kEval2_1temp * x + kEval2_2temp * x^2  
            }
            ##An monotonically increasing function will have diff(x) all > or equal to 0:
            if(all(diff(f2(Seq_Level)) >= 0)&length(LakeAreasel)>=10)
            {
                infl <- c(FALSE, diff(diff(f2(Lrange))>0)!=0)
                arr_infl<-Lrange[infl]
                lowlimittemp=lowlimit
                uplimittemp=uplimit
                if(length(arr_infl)==1)
                {
                   if(minL>arr_infl[1])
                   {
                      lowlimittemp=arr_infl[1]+0.11
                   }
                   else
                   {
                      uplimittemp=arr_infl[1]-0.11
                   }
                }
                if(length(Areaforvolume)>0)
                {
                    if((f2(lowlimittemp)<min(Areaforvolume))&(f2(uplimittemp)>max(Areaforvolume)))
                    {
                        Intercept = Intercept2_temp
                        kEval1 = kEval2_1temp
                        kEval2 = kEval2_2temp
                        f=f2
                        fit=fit2shift
                        lowlimit=lowlimittemp
                        uplimit=uplimittemp
                    }
                }
            }
        }
    }
    ##Try Cubic equation
    fit3shift<-lm(LakeAreasel~Level1+Level2+Level3)
    ##Using AIC to avoid overfitting
    if(length(coef(summary(fit3shift)))>12 & AIC(fit3shift)<AIC(fit) & !(isNslope)){
        kEval3_1temp <- coef(summary(fit3shift))["Level1", "Estimate"]    
        kEval3_2temp <- coef(summary(fit3shift))["Level2", "Estimate"] 
        kEval3_3temp <- coef(summary(fit3shift))["Level3", "Estimate"] 
        Intercept3_temp <- coef(summary(fit3shift))["(Intercept)", "Estimate"]
        Intercept3_temp=Intercept3_temp-kEval3_1temp*min(LakeLevelsel)-kEval3_2temp*min(LakeLevelsel^2)-kEval3_3temp*min(LakeLevelsel^3)
        if(!is.na(kEval3_3temp)&(summary(fit3shift)$r.squared>R2inHypsometry_thres))
        {
            f3 <- function(x) {
            Intercept3_temp + kEval3_1temp * x + kEval3_2temp * x^2 + kEval3_3temp * x^3
            }
            ##An monotonically increasing function will have diff(x) all > or equal to 0:
            if(all(diff(f3(Seq_Level)) >= 0)&length(LakeAreasel)>=10)
            {
                infl <- c(FALSE, diff(diff(f3(Lrange))>0)!=0)
                arr_infl<-Lrange[infl]
                lowlimittemp=-1000
                uplimittemp=8848
                if(length(arr_infl)==1)
                {
                   if(minL>arr_infl[1])
                   {
                      lowlimittemp=arr_infl[1]+0.11
                   }
                   else
                   {
                      uplimittemp=arr_infl[1]-0.11
                   }
                }
                if(length(arr_infl)==2)
                {
                    if(maxL<(arr_infl[1]-0.11))
                    {
                       uplimittemp=arr_infl[1]-0.11
                    }
                    else if(minL>(arr_infl[2]+0.11))
                    {
                       lowlimittemp=arr_infl[2]+0.11
                    }
                    else
                    {
                       uplimittemp=arr_infl[2]-0.11
                       lowlimittemp=arr_infl[1]+0.11
                    }
                }
                if(length(Areaforvolume)>0)
                {
                    if((f3(lowlimittemp)<min(Areaforvolume))&(f3(uplimittemp)>max(Areaforvolume)))
                    {
                        Intercept = Intercept3_temp
                        kEval1 = kEval3_1temp
                        kEval2 = kEval3_2temp
                        kEval3 = kEval3_3temp
                        f=f3
                        fit=fit3shift
                        uplimit=uplimittemp
                        lowlimit=lowlimittemp
                    }
                }
            }
        }
    }

    r_square<-summary(fit)$r.squared
    adjr_square=summary(fit)$adj.r.squared
    if(isNslope){
       r_square=0
       adjr_square=0
    }

    subname=paste(subname,round(r_square,3),sep=" (R^2:")
    subname=paste(subname,")",sep="")
    
    mydata.sel <- do.call(rbind, Map(data.frame, LakeLevelsel=LakeLevelsel, LakeAreasel=LakeAreasel))
    Hypsometry.fitted <- do.call(rbind, Map(data.frame, Seq_Level=Seq_Level, Areafitted=f(Seq_Level)))
    ##Plot Calibrated hypsometry
    p3 <- ggplot(mydata.sel, mapping = aes(x=LakeLevelsel, y=LakeAreasel)) + 
        ggtitle("Calibrated hypsometry") +
        xlab("Level (m)") +
        ylab("Water area (m2)") +
        geom_point() +
        geom_line(data = Hypsometry.fitted, aes(x=Seq_Level, y=Areafitted)) +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
        
    grid.arrange(p1, p2, p3, ncol = 1)
    
    summary(fit)$r.squared
    mse <- mean(residuals(fit)^2)
    rmse <- sqrt(mse)
    levelerrorOri = mydata$LevelErr
    levelerror=levelerrorOri[!is.na(levelerrorOri)]
    funSq <- function(x) {
        x * x
    }
    LerrorSq = funSq(levelerror)
    LerrorMean = sqrt(sum(LerrorSq))
    LerrorMean = LerrorMean/length(LerrorSq)
    
    ##Get Inverse of the function f
    inverse = function(f, lower = -100, upper = 100) {
        function(y) uniroot((function(x) f(x) - y), lower = lower, 
            upper = upper)[1]
    }
    f_inverse = inverse(f, lowlimit,uplimit)
    ##Output lake volume time series
    sink(outputfilefull)
    cat("Time,Mean_Levels,Mean_Level_Errors,LakeArea,NonCCover,Cleaned_Area,wsc,rws,,,")
    cat("Area=")
    cat(Intercept)
    if(kEval1>0){ cat("+") }
    cat(kEval1)
    cat("*")
    cat("x")
    if(abs(kEval2)>1e-8)
    {
       if(kEval2>0){ cat("+") }
       cat(kEval2)
       cat("*")
       cat("x^2")
    }
    if(abs(kEval3)>1e-8)
    {
       if(kEval3>0){ cat("+") }
       cat(kEval3)
       cat("*")
       cat("x^3")
    }

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
                cat(AreaOri[i]/Scalekm2_m2)
                cat(",")
                cat(NonCCover[i])
                if(!is.na(Area[i])){
                    cat(",")
                    cat(Area[i]/Scalekm2_m2)
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
                if (lastlevel < Level[i]) 
                  tempV = integrate(f, lastlevel, Level[i])$val else tempV = -1 * (integrate(f, Level[i], lastlevel)$val)
                cat(tempV/1e9)
                cat(",")
                cat(tempV/1e9 + lastV)#cat(tempV/1000 + lastV)
                lastV = tempV/1e9 + lastV #tempV/1000 + lastV
            } else {
                cat(",0,0")
                isFirstvalid = 1
            }
            lastlevel = Level[i]
        } else if (!is.na(Area[i])) {
            cat(f_inverse(Area[i])$root)
            cat(",,")
            cat(AreaOri[i]/Scalekm2_m2)
            cat(",")
            cat(NonCCover[i])
            cat(",")
            cat(Area[i]/Scalekm2_m2)
            if (i > 1 & isFirstvalid > 0) {
                cat(",")
                if (lastlevel < f_inverse(Area[i])$root) 
                  tempV = integrate(f, lastlevel, f_inverse(Area[i])$root)$val else tempV = -1 * (integrate(f, f_inverse(Area[i])$root, 
                  lastlevel)$val)
                cat(tempV/1e9)
                cat(",")
                cat(tempV/1e9 + lastV)
                lastV = tempV/1e9 + lastV
            } else {
                cat(",0,0")
                isFirstvalid = 1
            }
            lastlevel = f_inverse(Area[i])$root
        } else if (!is.na(AreaOri[i])) {
            cat(",,")
            cat(AreaOri[i]/Scalekm2_m2)
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

