###################################################################################
#####################Validating estimated lake water storage#######################
###################################################################################

library(mblm)
require(ggplot2)
library(mosaic)
library(stringr)
library(Metrics)
options("scipen"=100, "digits"=8)

##input: folder that contains in situ measurements (formatted)
MonthlyFormattedStorageDir <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/All/"
##input: folder that contains estimated lake water storage from satellites
MonthlyEstimateLWSDir <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Output/0Monthly/Volume_TS/All_1992_18/GlakesAllTypes/WaterVolume/0Reservoirs/"
##output folder for pairs of validation (csv file or a figure)
OutDir <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/"

##period for comparison
StartYear = 1992
EndYear = 2020

##Is plot only? If pairs of estimated and observated values being processed, may choose this option
IsPlotOnly=TRUE

########################Reservoir sedimentation adjustment if not applied yet#################################################################
##Reservoir info csv for reading the reservoir storage information
Reservoirinfocsv <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/reservoirs_LakeID_Capkm3.csv"
##Mean reservoir sedimentation rate
meanannuallossrate_sediment=-0.00245;
ReservoirInfo <- read.csv(file=Reservoirinfocsv)

cat("Cautions: problematic areas need to be updated and incoporated into the validation 211031undo\n")
Excluded_lakeidlist <- c("933")#,"678","738","1613","1068","629","1906","1335","2628","1008","3831","2174","2826")


if(!IsPlotOnly){
  filelist_insituVolume <- list.files(MonthlyFormattedStorageDir,'.csv',full.name = T)
  LakeIDlist_insitu <- as.numeric(gsub(".*?([0-9]+).*", "\\1", basename(filelist_insituVolume)))
  
  filelist_estimateLWS <- list.files(MonthlyEstimateLWSDir,'.csv',full.name = T)
  
  USArmyC <- list.files("C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/USArmy/",'.csv',full.name = T)
  LakeIDlist_USArmyC <- as.numeric(gsub(".*?([0-9]+).*", "\\1", basename(USArmyC)))
  
  data=matrix(scan("C:/0UCBShare/0Desktop/Papers/GlakesVolume/Data/ReservoirSedimentRates/US_reservoir_sedimentation_AnnualrateMatrix.txt"),ncol=2,byrow=T)
  LakeIDlist_USArmyC.sediment=data[,1]
  
  
  ####Set df for all pairs
  #Allpairs_df <- data.frame(matrix(ncol = 4, nrow = 0))
  #colnames <- c("Year", "Month", "Estimate","Insitu")
  #colnames(Allpairs_df) <- colnames
  Allpairs_df <- data.frame()
  Allpairs_df.sumbylake <- data.frame()
  for(estimateLWSfile in filelist_estimateLWS){
    lakeid <- as.numeric(gsub(".*?([0-9]+).*", "\\1", basename(estimateLWSfile)))
    if(!(lakeid %in% LakeIDlist_insitu) | lakeid %in% Excluded_lakeidlist)
    {
      next
    }
    LWSestimatefile <- read.csv(file=estimateLWSfile)
    LWSestimatefile$Year <- as.integer(substr(LWSestimatefile$Time,0,4))
    LWSestimatefile$Month <- as.integer((as.integer(substr(LWSestimatefile$Time,5,7))+20)/30+1)
    
    ReservoirInfoFind = ReservoirInfo[ReservoirInfo$LakeID==lakeid,]
    annualloss.sedimentation = 0
    
    if(dim(ReservoirInfoFind)[1]>0 & ReservoirInfoFind$NReservoir[1]==0)
      annualloss.sedimentation = meanannuallossrate_sediment*ReservoirInfoFind$Cap_km3[1]
    
    
    
    insitufile <- paste0(MonthlyFormattedStorageDir,lakeid,".csv")
    insitu_data <- read.csv(file=insitufile)
    insitu_data_filtered <- insitu_data[insitu_data$YEAR>=StartYear & insitu_data$YEAR<=EndYear, ]
    UniqueYears <- unique(insitu_data_filtered$YEAR)
    insitu_data_na <- insitu_data[is.na(insitu_data$YEAR),]
    if(dim(insitu_data_na)[1]>0)
    {
      cat(paste0(basename(insitufile)," year column has NA value; please check","\n"))
      next
    }
    if((max(UniqueYears)-min(UniqueYears)+1)>length(UniqueYears))
    {
      cat(paste0(basename(insitufile)," missing values in at least one year; please check","\n"))
      next
    }
    data_mergetemp <- merge(insitu_data_filtered, LWSestimatefile, by.x = c('YEAR','MONTH'), by.y = c('Year','Month'))
    data_mergetemp_filter <- data_mergetemp[!is.na(data_mergetemp$meanLWS) & !is.na(data_mergetemp$rws),]
    
    ###########Adjustment for sedimenation
    data_mergetemp_filter$rws = data_mergetemp_filter$rws + (data_mergetemp_filter$YEAR-StartYear)*annualloss.sedimentation
    
    if(dim(data_mergetemp_filter)[1]==0){
      next
    }
    data_mergetemp_filter$meanLWS <- data_mergetemp_filter$meanLWS - mean(data_mergetemp_filter$meanLWS)
    data_mergetemp_filter$rws <- data_mergetemp_filter$rws - mean(data_mergetemp_filter$rws)
    Allpairs_df <- rbind(Allpairs_df, data.frame(LakeID=rep(lakeid,length(data_mergetemp_filter$YEAR)),Year=data_mergetemp_filter$YEAR,Month=data_mergetemp_filter$MONTH,
                                                 Estimate=data_mergetemp_filter$rws,Insitu=data_mergetemp_filter$meanLWS))
    cat(lakeid)
    cat("\n")
    #plot(data_mergetemp_filter$meanLWS,data_mergetemp_filter$rws)
    Allpairs_df.sumbylake <- rbind(Allpairs_df.sumbylake, cbind.data.frame(LakeID = lakeid,
                                                                           Number = length(data_mergetemp_filter$meanLWS),
                                                                           RMSE = rmse(data_mergetemp_filter$meanLWS,data_mergetemp_filter$rws)))
  }
  write.csv(Allpairs_df,paste(OutDir,"Glakes_Validation_pairs.csv"), row.names = F)
  write.csv(Allpairs_df.sumbylake,paste(OutDir,"Glakes_Validation_summarybyeachlake.csv"), row.names = F)
}else{
  Allpairs_df=read.csv(file=paste(OutDir,"Glakes_Validation_pairs.csv"), header=TRUE, sep=",")
}

##Get statistic metrics round(corcheck$estimate[1],3)
cat("Validate estimate LWS anomalies against in situ observations\n")
fit<- lm(Allpairs_df$Insitu~Allpairs_df$Estimate) 
cat(paste0("A total of ",length(Allpairs_df$Insitu)," independent in situ observations for comparison\n"))
All.r_square<-summary(fit)$r.squared
cat(paste0("R-square:",round(All.r_square,2),"\n"))
All.RMSE = rmse(Allpairs_df$Insitu,Allpairs_df$Estimate)
cat(paste0("RMSE:",round(All.RMSE,2),"\n"))


jpeg(paste(OutDir,"Validation_LWS_estimated_vs_insitumeasurements.jpg"), width = 8, height = 8, units = 'in', res = 600) 
ggplot(Allpairs_df)+
  geom_point(aes(x = Insitu, y = Estimate, color = 'red'), size = 3, alpha = 1, show.legend = FALSE)+
  xlim(c(-10,13))+
  ylim(c(-10, 13))+
  xlab("In situ storage anomaly (Gt)")+
  ylab("Estimated storage anomaly (Gt)")+
  theme_bw()+
  theme(text = element_text(size = 24))+
  theme(legend.position =  c(0.2,0.9),
        legend.title = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.text = element_text(color = 'black'))
dev.off()
