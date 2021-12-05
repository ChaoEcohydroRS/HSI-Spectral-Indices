############################################################################################################
##Constructing multiple regression ensemble to diagnose the dominant driver of lake volume changes##########
##consider hydroclimate and human variables and model lake volume changes over 1992-2010####################
############################################################################################################
setwd("G:/My Drive/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Glakes/")

library(stringr)
library(stringr)
library(MASS)
library(tools)

##input: folder that contains yearly lake volume stored in a csv format
LakeVolumeDir="./data/Output_Volume_TimeSeries/Yearly/"
##input: root folder for data of hydroclimate and human variables
ClimateDataRootFolder="./data/Yearly_Drivers_TimeSeries/"
##outputfolder: multiple regression ensemble for lake drivers
OutputFolder= "./data/Output_Modelensembel/"

cat("Creat subfolders in the ClimateDataRootFolder and input datasets of yearly Prep, PET, Tmean, Runoff, HW as csv files in each subfolder\n")
Prep.DS <- c("PtrMSWEP","PtrERA","PtrMERRA","PtrCRU")
PET.DS <- c("PEtrGLEAM","PEtrCRU")
Tmean.DS <- c("TtrCRU","TtrERA","TtrGHCN")
##Human water consumption data based on four hydrologic models and biased corrected (Huang et al. 2018, HESS)
WC.DS <- c("WCtrH08_","WCtrJPJ","WCtrPCR","WCtrWGAP")
RunoffFoldername="RtrGRADES"

##Set the period for modeling lake volume change
#cat("For north hemisphere using 1993 (199210-199309); for south hemisphere using 1994(199304-199403)")
#For Oct: use 1993; For Apr: use 1994
StartYear = 1992
EndYear =2010

##whether scale variables in the regression model
isScaleData=TRUE

##Adding details of the info in output file names
outfilerootname_DetailsLM="0wWateruse_BIC_Details_StepwiseLM_"
outfilerootname_SummaryLM="0wWateruse_BIC_Sum_StepwiseLM_"
outfile_datasourcecomparison="0wWateruse_BIC_Sum_ModelPerformancebyDataSource"

################Columns of csv files to be read####################
##Year column
YearFieldinAllcsv="Year"
##Lake volume change column
ValueFieldinLakevolume="VolumeChange"
##Column for values of each explanatory variable
ValueFieldinAllcsvexpectlakevolume="Value"

##To record whether a lake has a large volume rate 
LargeVolumerateThres = 0.1 ##Gt per year

##Switch to the 2nd year of the study period in order to calculate lake volume change 
StartYear = StartYear +1

##func for extracting the first number from a string
numextract <- function(string){ 
    str_extract(string, "\\-*\\d+\\.*\\d*")
} 

cat("Premininary examination of the relationship between volume change and each explanatory variable\n")
cat("Example for Aral Sea: links to human water use\n")
LakeIDToBeExamined="82"
DriverToBeExamined=WC.DS
##Plot setting 4 panels: 2 rows * 2 columns
par(mfrow = c(2, 2))
par(mar=c(2, 2, 2, 2))
for(eachdriverfoldercheck in DriverToBeExamined)
{
    inputfilecheck <- paste0("ID",LakeIDToBeExamined,".csv")
    LakeVolumefilecheck = paste0(LakeVolumeDir, file_path_sans_ext(inputfilecheck), "_yearly.csv")
    Lakevolumedatacheck <- read.csv(file=LakeVolumefilecheck, header=TRUE, sep=",")
    LWScheck=Lakevolumedatacheck[Lakevolumedatacheck[YearFieldinAllcsv]>=StartYear & Lakevolumedatacheck[YearFieldinAllcsv]<=EndYear,]$VolumeChange
    
    DriverDircheck = paste0(ClimateDataRootFolder,eachdriverfoldercheck,"/")
    Driverfilecheck = paste(DriverDircheck, inputfilecheck, sep = "")
    Driverdatacheck <- read.csv(file=Driverfilecheck, header=TRUE, sep=",")
    
    Drivercheck=Driverdatacheck[Driverdatacheck[YearFieldinAllcsv]>=StartYear & Driverdatacheck[YearFieldinAllcsv]<=EndYear,]$Value
    
    fitcheck<- lm(LWScheck~Drivercheck) 
    kEvalcheck <- coef(summary(fitcheck))["Drivercheck", "Estimate"]
    Interceptcheck <- coef(summary(fitcheck))["(Intercept)", "Estimate"]
    fcheck <- function(x) {
        Interceptcheck + kEvalcheck * x 
    }
    corcheck <-cor.test(Drivercheck, LWScheck,method = "pearson")
    titlecheck=paste0(LakeIDToBeExamined," Cor:",round(corcheck$estimate[1],3),eachdriverfoldercheck)
    plot(Drivercheck,LWScheck,main=titlecheck,xlab="driver variable",ylab="Lake volume change")
    Seq_Levelcheck <- seq(from=min(Drivercheck), to=max(Drivercheck), by=0.01)
    lines(Seq_Levelcheck,fcheck(Seq_Levelcheck))
}

##construct regression models from different data sources
func_lakedriver <- function(driversDataSourcestr){
    comarr <- unlist(strsplit(driversDataSourcestr, "[|]"))
    PrepDir=paste0(ClimateDataRootFolder,comarr[1],"/")
    PETDir=paste0(ClimateDataRootFolder,comarr[2],"/")
    RiverrunoffDir=paste0(ClimateDataRootFolder,RunoffFoldername,"/")
    TmeanDir=paste0(ClimateDataRootFolder,comarr[3],"/")
    WCDir=paste0(ClimateDataRootFolder,comarr[4],"/")
    ##set output to be written
    outputdetailstxt=paste0(OutputFolder,outfilerootname_DetailsLM,str_replace_all(com,"[|]","_"),".txt")
    outputsummarycsv=paste0(OutputFolder,outfilerootname_SummaryLM,str_replace_all(com,"[|]","_"),".csv")
    fileList <- list.files(PrepDir,'.csv$', full.name = F)
    sink(outputdetailstxt)
    cat("LakeID Volumerate Volumeratepval LMRsquared FStatistic FStatisticpval |ModelFormula|ModelCoefficients|ModelCoefficientspval|MeanSquaresbyeachvariable\n")
    vec_lakeid <- c()
    vec_srate <- c()
    vec_rsquaredAll <- c()
    vec_rsquaredLarge <- c()
    vec_fstatisticAll <- c()
    vec_fstatisticLarge <- c()
    vec_fstatisticAllpval <- c()
    vec_fstatisticLargepval <- c()
    for(inputfile in fileList)
    {
        prepfile = paste(PrepDir, inputfile, sep = "")
        PETfile = paste(PETDir, inputfile, sep = "")
        Runofffile = paste(RiverrunoffDir, inputfile, sep = "")
        WCfile = paste(WCDir, inputfile, sep = "")
        Tmeanfile = paste(TmeanDir,inputfile,sep = "")
        LakeVolumefile = paste0(LakeVolumeDir, file_path_sans_ext(inputfile), "_yearly.csv")
        ##skip if one of the explainatory variables does not have data
        if(file.exists(PETfile)==FALSE|file.exists(Runofffile)==FALSE|
           file.exists(LakeVolumefile)==FALSE|file.exists(WCfile)==FALSE)
        {
            next
        }
        
        LakeID <- numextract(substr(file_path_sans_ext(inputfile),3, nchar(file_path_sans_ext(inputfile))))
        Lakevolumedata <- read.csv(file=LakeVolumefile, header=TRUE, sep=",")
        LWS=Lakevolumedata[Lakevolumedata[YearFieldinAllcsv]>=StartYear & Lakevolumedata[YearFieldinAllcsv]<=EndYear,]$VolumeChange
        ##Get lake volume trend
        Srate <- Lakevolumedata$VolumeRate[1]
        ##Get pvalue of the trend 
        Spvalue=Lakevolumedata$VolumePvalue[1]
        cat(LakeID)
        cat(" ")
        cat(Srate)
        cat(" ")
        cat(Spvalue)
        cat(" ")
        Prepdata <- read.csv(file=prepfile, header=TRUE, sep=",")
        Prep=Prepdata[Prepdata[YearFieldinAllcsv]>=StartYear & Prepdata[YearFieldinAllcsv]<=EndYear,]$Value
        PETdata <- read.csv(file=PETfile, header=TRUE, sep=",")
        PET=PETdata[PETdata[YearFieldinAllcsv]>=StartYear & PETdata[YearFieldinAllcsv]<=EndYear,]$Value
        Runoffdata <- read.csv(file=Runofffile, header=TRUE, sep=",")
        Runoff=Runoffdata[Runoffdata[YearFieldinAllcsv]>=StartYear & Runoffdata[YearFieldinAllcsv]<=EndYear,]$Value
        Tmeandata <- read.csv(file=Tmeanfile, header=TRUE, sep=",")
        Tmean <- Tmeandata[Tmeandata[YearFieldinAllcsv]>=StartYear & Tmeandata[YearFieldinAllcsv]<=EndYear,]$Value
        WCdata <- read.csv(file=WCfile, header=TRUE, sep=",")
        WC=WCdata[WCdata[YearFieldinAllcsv]>=StartYear & WCdata[YearFieldinAllcsv]<=EndYear,]$Value
        
        #ID1837Goud_eZareh dried out by 2003
        if(LakeID==1837)
        {
            Prep=Prepdata[Prepdata[YearFieldinAllcsv]>=StartYear & Prepdata[YearFieldinAllcsv]<=2002,]$Value
            PET=PETdata[PETdata[YearFieldinAllcsv]>=StartYear & PETdata[YearFieldinAllcsv]<=2002,]$Value
            Runoff=Runoffdata[Runoffdata[YearFieldinAllcsv]>=StartYear & Runoffdata[YearFieldinAllcsv]<=2002,]$Value
            Tmean <- Tmeandata[Tmeandata[YearFieldinAllcsv]>=StartYear & Tmeandata[YearFieldinAllcsv]<=2002,]$Value
            WC=WCdata[WCdata[YearFieldinAllcsv]>=StartYear & WCdata[YearFieldinAllcsv]<=2002,]$Value
            LWS=Lakevolumedata[Lakevolumedata[YearFieldinAllcsv]>=StartYear & Lakevolumedata[YearFieldinAllcsv]<=2002,]$VolumeChange
        }
        #ID7386Toshka filled by human act in Nov 1998 and dried out by 2015
        if(LakeID==7386)
        {
            Prep=Prepdata[Prepdata[YearFieldinAllcsv]>=2000 & Prepdata[YearFieldinAllcsv]<=EndYear,]$Value
            PET=PETdata[PETdata[YearFieldinAllcsv]>=2000 & PETdata[YearFieldinAllcsv]<=EndYear,]$Value
            Runoff=Runoffdata[Runoffdata[YearFieldinAllcsv]>=2000 & Runoffdata[YearFieldinAllcsv]<=EndYear,]$Value
            Tmean <- Tmeandata[Tmeandata[YearFieldinAllcsv]>=2000 & Tmeandata[YearFieldinAllcsv]<=EndYear,]$Value
            WC=WCdata[WCdata[YearFieldinAllcsv]>=2000 & WCdata[YearFieldinAllcsv]<=EndYear,]$Value
            LWS=Lakevolumedata[Lakevolumedata[YearFieldinAllcsv]>=2000 & Lakevolumedata[YearFieldinAllcsv]<=EndYear,]$VolumeChange
        }
        #ID8 Great Bear, volume since 200201
        if(LakeID==8)
        {
            Prep=Prepdata[Prepdata[YearFieldinAllcsv]>=2003 & Prepdata[YearFieldinAllcsv]<=EndYear,]$Value
            PET=PETdata[PETdata[YearFieldinAllcsv]>=2003 & PETdata[YearFieldinAllcsv]<=EndYear,]$Value
            Runoff=Runoffdata[Runoffdata[YearFieldinAllcsv]>=2003 & Runoffdata[YearFieldinAllcsv]<=EndYear,]$Value
            Tmean <- Tmeandata[Tmeandata[YearFieldinAllcsv]>=2003 & Tmeandata[YearFieldinAllcsv]<=EndYear,]$Value
            WC=WCdata[WCdata[YearFieldinAllcsv]>=2003 & WCdata[YearFieldinAllcsv]<=EndYear,]$Value
            LWS=Lakevolumedata[Lakevolumedata[YearFieldinAllcsv]>=2003 & Lakevolumedata[YearFieldinAllcsv]<=EndYear,]$VolumeChange
        }

        WCstd <- sd(WC, na.rm = FALSE)
        LWCstd <- sd(LWS, na.rm = FALSE)
        if(isScaleData)
        {
            LWS=scale(LWS)
            Prep=scale(Prep)
            PET=scale(PET)
            Runoff=scale(Runoff)
            Tmean=scale(Tmean)
            WC=scale(WC)
        }
        
        vec_indvar <- c()
        if(!is.na(Prep[1]) & max(Prep)!=min(Prep))
        {
            vec_indvar=append(vec_indvar,"Prep")
        }
        if(!is.na(PET[1]) & max(PET)!=min(PET))
        {
            vec_indvar=append(vec_indvar,"PET")
        }
        if(!is.na(Runoff[1]) & max(Runoff)!=min(Runoff))
        {
            vec_indvar=append(vec_indvar,"Runoff")
        }
        if(!is.na(Tmean[1]) & max(Tmean)!=min(Tmean))
        {
            vec_indvar=append(vec_indvar,"Tmean")
        }
        if(!is.na(WC[1]) & max(WC)!=min(WC) & WCstd > 0.05*LWCstd)
        {
            vec_indvar=append(vec_indvar,"WC")
        }
        if(length(vec_indvar)==0)
        {
            next
        }
        lmformulastr <- vec_indvar[1]
        Isskipfirst <- TRUE
        for(eachvar in vec_indvar)
        {
            if(Isskipfirst)
            {
                Isskipfirst<-FALSE
                next
            }
            lmformulastr=paste(lmformulastr,eachvar,sep="+")
        }
        lmformulastr<-paste("LWS",lmformulastr,sep="~")
        ####Using AIC
        ##stepAIC(lm(formula(lmformulastr)),direction = "both",trace = FALSE)
        ####Using BIC
        ##step(lm(formula(lmformulastr)),direction ="both",k=log(length(Prep)),trace = FALSE)
        ##Using BIC to avoid overfitting as BIC favors parsimony
        stepwisemodel <- step(lm(formula(lmformulastr)),direction ="both",k=log(length(Prep)),trace = FALSE)
        
        formulastr=toString(formula(stepwisemodel))
        poss=regexpr(",",formulastr)[1]
        formulastr=substr(formulastr,poss[1]+1,nchar(formulastr))
        formulastr=str_replace(formulastr,"LWS,","|")
        if(dim(coef(summary(stepwisemodel)))[1]<=1)
        {
            cat("NA NA | NA\n")
            next
        }
        R2 <- summary(stepwisemodel)$r.squared
        fstatistic <- summary(stepwisemodel)$fstatistic[1]
        fstatistic.pval <- pf(summary(stepwisemodel)$fstatistic[1],summary(stepwisemodel)$fstatistic[2],summary(stepwisemodel)$fstatistic[3],lower.tail=FALSE)
        cat(R2)
        cat(" ")
        cat(fstatistic)
        cat(" ")
        cat(fstatistic.pval)
        cat(" ")
        cat(formulastr)
        fit<-lm(formula(stepwisemodel))
        cat("|")
        cat(coef(summary(fit))[, "Estimate"][2:length(coef(summary(fit))[, "Estimate"])])
        cat("|")
        cat(coef(summary(fit))[, "Pr(>|t|)"][2:length(coef(summary(fit))[, "Pr(>|t|)"])])
        cat("|")
        anova(fit)
        cat(anova(fit)$"Mean Sq"/sum(anova(fit)$"Mean Sq"))
        cat("\n")
        ##for writing a log if needed
        if(abs(Srate)>0.1)
        {
            vec_rsquaredLarge=append(vec_rsquaredLarge,R2)
            vec_fstatisticLarge=append(vec_fstatisticLarge,fstatistic)
            vec_fstatisticLargepval=append(vec_fstatisticLargepval,fstatistic.pval)
        }
        vec_lakeid=append(vec_lakeid,LakeID)
        vec_srate=append(vec_srate,Srate)
        vec_rsquaredAll=append(vec_rsquaredAll,R2)
        vec_fstatisticAll=append(vec_fstatisticAll,fstatistic)
        vec_fstatisticAllpval=append(vec_fstatisticAllpval,fstatistic.pval)
    }
    sink()
    #sink(outputsummarycsv)
    #cat("LakeID,Rate,Rsquared,Fstatistics,Fstatisticspval,IsLargeRate,\n")
    #for(bgi in 1:length(vec_rsquaredAll))
    #{
    #     cat(vec_lakeid[bgi])
    #     cat(",")
    #     cat(vec_srate[bgi])
    #     cat(",")
    #     cat(vec_rsquaredAll[bgi])
    #     cat(",")
    #     cat(vec_fstatisticAll[bgi])
    #     cat(",")
    #     cat(vec_fstatisticAllpval[bgi])
    #     cat(",")
    #     cat(as.numeric(abs(vec_srate[bgi])>LargeVolumerateThres))
    #     cat(",\n")
    # }
    # sink()
    return(paste0(length(vec_rsquaredAll),",",mean(vec_rsquaredAll),",",mean(vec_fstatisticAll),",",
                  length(vec_rsquaredLarge),",",mean(vec_rsquaredLarge),",",mean(vec_fstatisticLarge),","))
}

##Initialize for manually checking when needed 
Prepdata <- Prep.DS[1]
PETdata <- PET.DS[1]
Tmeandata <- Tmean.DS[1]
WCdata <- WC.DS[1]

vec_statisticsummarystr <- c()
for(Prepdata in Prep.DS)
{
    for(PETdata in PET.DS)
    {
            for(WCdata in WC.DS)
            {
                for(Tmeandata in Tmean.DS)
                {
                    com <- paste0(Prepdata,"|",PETdata,"|",Tmeandata,"|",WCdata)
                    statisticsummarystr <- func_lakedriver(com)
                    vec_statisticsummarystr=append(vec_statisticsummarystr,paste0(statisticsummarystr,com,","))
                }
            }
    }
}
cat("All models were run completely!Please check model outputs\n")
#sink(paste0(ClimateDataRootFolder,outfile_datasourcecomparison,".csv"))
#cat("NumAll,RsquaredAll,FstatisticAll,NumLargerate,RsquaredLargerate,FstatisticLargerate,\n")
#for(statstr in vec_statisticsummarystr)
#{
#    cat(statstr)
#    cat("\n")
#}
#sink()


