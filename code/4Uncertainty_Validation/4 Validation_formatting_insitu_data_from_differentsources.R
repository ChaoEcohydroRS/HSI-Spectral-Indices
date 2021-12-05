###################################################################################
###########Formatting In situ storage data from different sources #################
###################################################################################
library(mblm)
require(ggplot2)
library(mosaic)
library(stringr)
options("scipen"=100, "digits"=8)

###########Insitu U.S. Army Corps Begin#################
InsituStorageDir <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/USArmyCorps/Corps Data Shared/Measurements/wStorage/"
Outputfolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/USArmy/0/"
#Data from U.S. Army Corps use NID as identifier while our lake dataset uses LakeID as identifier
MatchNID_LakeIDfile <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/Match_NIDID_LakeID.csv"
NID_LakeIDpair_data <- read.csv(file=MatchNID_LakeIDfile)
NIDID_filtered <- NID_LakeIDpair_data$NIDID

filelist_insituVolume <- list.files(InsituStorageDir,'.csv',full.name = T)

for(insitufile in filelist_insituVolume){
    insitu_data <- read.csv(file=insitufile)
    insitu_data_filtered <- insitu_data[insitu_data$NIDID %in% NIDID_filtered & !is.na(insitu_data$STORAGE_ACFT), ]
    #if id is different by= c("ID1","ID2")
    insitu_data_filtered <- merge(insitu_data_filtered,NID_LakeIDpair_data, by="NIDID")
    insitu_data_filteredyearly <- insitu_data_filtered %>%
           group_by(LakeID,YEAR,MONTH) %>%
           summarise(meanLWS = mean(STORAGE_ACFT, na.rm=T)*0.000001233482)
    UniqueLakeIDs <- unique(insitu_data_filteredyearly$LakeID)
    for(eachID in UniqueLakeIDs){
        outcsv <- paste0(Outputfolder,eachID,".csv")
        write.csv(insitu_data_filteredyearly[insitu_data_filteredyearly$LakeID==eachID,][c("YEAR","MONTH","meanLWS")], outcsv, row.names = F)
    }
}
###########Insitu U.S. Army Corps End#################

###########Insitu Australia Begin#################
Datafolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/Australia/"
Outputfolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/AU/"
filelist_insituVolume <- list.files(Datafolder,'.csv',full.name = T)
for(insitufile in filelist_insituVolume){
    insitu_data_filtered <- read.csv(file=insitufile,header = FALSE, comment.char = "#")
    insitu_data_filtered$YEAR = as.integer(substr(insitu_data_filtered$V1,1,4))
    insitu_data_filtered$MONTH= as.integer(substr(insitu_data_filtered$V1,6,7))
    insitu_data_filteredyearly <- insitu_data_filtered[!is.na(insitu_data_filtered$V2),] %>%
        group_by(YEAR,MONTH) %>%
        summarise(meanLWS = mean(V2, na.rm=T)*0.000001)
    
    outcsv <- paste0(Outputfolder,as.numeric(gsub(".*?([0-9]+).*", "\\1", basename(insitufile))),".csv" )
    write.csv(insitu_data_filteredyearly, outcsv, row.names = F)
}
###########Insitu Australia End#################

###########Insitu Texas Water Begin#################
Datafolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/TexasWater/"
Outputfolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/TexasWater/"
filelist_insituVolume <- list.files(Datafolder,'ID',full.name = T)
for(insitufile in filelist_insituVolume){
    insitu_data_filtered <- read.csv(file=insitufile,header = TRUE, comment.char = "#")
    insitu_data_filtered$YEAR = as.integer(substr(insitu_data_filtered$date,1,4))
    insitu_data_filtered$MONTH= as.integer(substr(insitu_data_filtered$date,6,7))
    insitu_data_filteredyearly <- insitu_data_filtered[!is.na(insitu_data_filtered$reservoir_storage),] %>%
        group_by(YEAR,MONTH) %>%
        summarise(meanLWS = mean(reservoir_storage, na.rm=T)*0.000001233482)
    outcsv <- paste0(Outputfolder,as.numeric(gsub(".*?([0-9]+).*", "\\1", basename(insitufile))),".csv" )
    write.csv(insitu_data_filteredyearly, outcsv, row.names = F)
}
###########Insitu Texas Water End#################

###########Insitu California Water Begin#################
Datafolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/CaliforniaWater/rawdata_org/"
Outputfolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/CaliforniaWater/"
MatchCAWaterID_LakeIDfile <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/Match_CAWaterID_LakeID.csv"
CAWaterID_LakeIDpair_data <- read.csv(file=MatchCAWaterID_LakeIDfile)
CAWaterID_filtered <- CAWaterID_LakeIDpair_data$ID

filelist_insituVolume <- list.files(Datafolder,'.csv',full.name = T)
for(insitufile in filelist_insituVolume){
    pos_datetemp <- str_locate(basename(insitufile),'_')[1]+1
    CAWaterIDtemp <- substr(basename(insitufile), pos_datetemp, pos_datetemp+2)
    if(!(CAWaterIDtemp %in% CAWaterID_filtered))
    {
        next
    }
    LakeIDtemp <- CAWaterID_LakeIDpair_data[CAWaterID_LakeIDpair_data$ID == CAWaterIDtemp,]$LakeID[1]
    insitu_data <- read.csv(file=insitufile)
    insitu_data_filtered <- insitu_data[!insitu_data$Value %in% c("---") & !is.na(insitu_data$Value), ]
    insitu_data_filtered$YEAR = as.integer(substr(insitu_data_filtered$OBS_Date,1,4))
    insitu_data_filtered$MONTH= as.integer(substr(insitu_data_filtered$OBS_Date,5,6))
    insitu_data_filtered$Value=as.double(insitu_data_filtered$Value)
    newdf <- data.frame(YEAR=insitu_data_filtered$YEAR, MONTH=insitu_data_filtered$MONTH,meanLWS=insitu_data_filtered$Value*0.000001233482)
    outcsv <- paste0(Outputfolder,LakeIDtemp,".csv")
    write.csv(newdf, outcsv, row.names = F)
}
###########Insitu California Water End#################


###########Insitu USGS Begin#################
Datafolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/USGS/ReservoirStorage/"
Outputfolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/USGS/"
filelist_insituVolume <- list.files(Datafolder,'ID',full.name = T)
for(insitufile in filelist_insituVolume){
  pos_datetemp <- str_locate(basename(insitufile),'_')[1]-1
  
  insitu_data_filtered <- read.csv(file=insitufile,header = TRUE, sep ="\t", comment.char = "#")
  insitu_data_filtered$YEAR = insitu_data_filtered$year_nu
  insitu_data_filtered$MONTH= insitu_data_filtered$month_nu 

  insitu_data_filteredmonthly <- insitu_data_filtered[!is.na(insitu_data_filtered$mean_va),] %>%
    group_by(YEAR,MONTH) %>%
    summarise(meanLWS = mean(mean_va, na.rm=T)*0.000001233482)
  
  outcsv <- paste0(Outputfolder,substr(basename(insitufile),3,pos_datetemp),".csv")
  write.csv(insitu_data_filteredmonthly, outcsv, row.names = F)
}
###########Insitu USGS End#################

###########Insitu Spain Begin#################
Datafolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/Spain/"
Outputfolder <- "C:/0D/Research/3GLakeWatch/Glakes_8dVolume/Validation/0Processed/MonthlyLWS/Spain/"
filelist_insituVolume <- list.files(Datafolder,'ID',full.name = T)
for(insitufile in filelist_insituVolume){
  pos_datetemp <- str_locate(basename(insitufile),'_')[1]-1
  
  insitu_data_filtered <- read.csv(file=insitufile,header = FALSE, sep =",")
  insitu_data_filtered$Storage <- as.double(str_replace_all(as.character(insitu_data_filtered$V2),",","")) 
  insitu_data_filtered <- insitu_data_filtered[!is.na(insitu_data_filtered$Storage),]
  insitu_data_filtered$YEAR = as.integer(substr(insitu_data_filtered$V1,nchar(insitu_data_filtered$V1)-6,nchar(insitu_data_filtered$V1)-2))
  insitu_data_filtered$MONTH= as.integer(substr(insitu_data_filtered$V1,nchar(insitu_data_filtered$V1)-1,nchar(insitu_data_filtered$V1)))

  insitu_data_filteredmonthly <- insitu_data_filtered[!is.na(insitu_data_filtered$Storage),] %>%
    group_by(YEAR,MONTH) %>%
    summarise(meanLWS = mean(Storage, na.rm=T)*0.001)
  
  outcsv <- paste0(Outputfolder,substr(basename(insitufile),3,pos_datetemp),".csv")
  write.csv(insitu_data_filteredmonthly, outcsv, row.names = F)
}
###########Insitu Spain End#################
