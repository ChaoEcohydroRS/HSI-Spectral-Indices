##############################################################################################
##Estimating the total population depending on lakes with a significant water loss (p<0.05)###
##############################################################################################

library(ncdf4)
library(Thermimage)
library(raster)
#library(doParallel)
#library(foreach)
library(lubridate)
#library(StationSWERegressionV2)
library("readxl")
library(mosaic)
#library(pryr)
library(hrbrthemes)
library(viridis)

##########################################Step1: Get Used Hydro ID (COMID)#################
Hydrolakeidcsv_fromlakeshp <- "/Users/fangfangyao/Documents/Proj/GlakesVolume/GLCP/Input/HydroLakes_GT100_Natural_Point_Q_P_PET_T_North_A_South_Hemi_wReservoirs.csv"
#lakeid for dryland lakes
lakeidfile_arid <- "/Users/fangfangyao/Documents/Proj/GlakesVolume/GLCP/Input/Hylak_id_drylands.csv"
#GLobal lake area, climate and population dataset
GLCPfile <- "/Users/fangfangyao/Documents/Proj/GlakesVolume/GLCP/Input/glcp.csv"
Outfigure <- "/Users/fangfangyao/Documents/Proj/GlakesVolume/GLCP/Output/Figure5_AffectedPopulation_wEnergy_Switch.jpg"
#Which year to use for counting population
Year4population <- 2015
nc_initialrunflg=FALSE

intermediateGLCPinYear <- paste0(substr(GLCPfile, 1, nchar(GLCPfile)-4),'_',Year4population,'.csv')

if(nc_initialrunflg){
  GLCPdata <- read.csv(file=GLCPfile)
  GLCPdata_filtered <- GLCPdata[GLCPdata$year ==Year4population, ]
  write.csv(GLCPdata_filtered, intermediateGLCPinYear)
}else
{
  GLCPdata_filtered=read.csv(intermediateGLCPinYear,sep = ',')
}

Hydrolakeid_candidate <- read.csv(Hydrolakeidcsv_fromlakeshp)

Hydrolakeid_freshwater <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID!=1,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_globalfreshwater  <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_freshwater$Hylak_id, ]
Hydrolakeid_envidegradation <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID!=3,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_envidegradation <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_envidegradation$Hylak_id, ]
Hydrolakeid_energyreduction <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID==3 & Hydrolakeid_candidate$Hydropower==1,]
GLCPdata_filtered_energyreduction <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_energyreduction$Hylak_id, ]
Hydrolakeid_globalAll <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_globalAll  <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_globalAll$Hylak_id, ]
#Hydrolakeid_bothimpacts <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID==0,]#& Hydrolakeid_candidate$TypeID==1
#GLCPdata_filtered_bothimpacts <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_bothimpacts$Hylak_id, ]

sum(GLCPdata_filtered_globalfreshwater$pop_sum)
sum(GLCPdata_filtered_envidegradation$pop_sum)
sum(GLCPdata_filtered_energyreduction$pop_sum)
sum(GLCPdata_filtered_globalAll$pop_sum)

Hydrolakeid_candidate <- read.csv(Hydrolakeidcsv_fromlakeshp)
lakeid_arid <- read.csv(lakeidfile_arid)
lakeid_arid <- lakeid_arid[lakeid_arid$Hylak_id>0,]
Hydrolakeid_candidate <- Hydrolakeid_candidate[Hydrolakeid_candidate$Hylak_id %in% lakeid_arid$Hylak_id, ]
Hydrolakeid_drylandfreshwater <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID!=1,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_drylandfreshwater  <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_drylandfreshwater$Hylak_id, ]
Hydrolakeid_drylandenvidegradation <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID!=3,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_drylandenvidegradation <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_drylandenvidegradation$Hylak_id, ]
Hydrolakeid_drylandenergyreduction <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID==3 & Hydrolakeid_candidate$Hydropower==1,]
GLCPdata_filtered_drylandenergyreduction <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_drylandenergyreduction$Hylak_id, ]
Hydrolakeid_drylandAll <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_drylandAll  <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_drylandAll$Hylak_id, ]
#Hydrolakeid_drylandbothimpacts <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID==0,]#& Hydrolakeid_candidate$TypeID==1
#GLCPdata_filtered_drylandbothimpacts <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_drylandbothimpacts$Hylak_id, ]

sum(GLCPdata_filtered_drylandfreshwater$pop_sum)
sum(GLCPdata_filtered_drylandenvidegradation$pop_sum)
sum(GLCPdata_filtered_drylandenergyreduction$pop_sum)
sum(GLCPdata_filtered_drylandAll$pop_sum)

Hydrolakeid_candidate <- read.csv(Hydrolakeidcsv_fromlakeshp)
lakeid_arid <- read.csv(lakeidfile_arid)
lakeid_arid <- lakeid_arid[lakeid_arid$Hylak_id>0,]
Hydrolakeid_candidate <- Hydrolakeid_candidate[!Hydrolakeid_candidate$Hylak_id %in% lakeid_arid$Hylak_id, ]
Hydrolakeid_wetterlandfreshwater <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID!=1,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_wetterlandfreshwater  <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_wetterlandfreshwater$Hylak_id, ]
Hydrolakeid_wetterlandenvidegradation <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID!=3,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_wetterlandenvidegradation <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_wetterlandenvidegradation$Hylak_id, ]
Hydrolakeid_wetterlandenergyreduction <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID==3 & Hydrolakeid_candidate$Hydropower==1,]
GLCPdata_filtered_wetterlandenergyreduction <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_wetterlandenergyreduction$Hylak_id, ]
Hydrolakeid_wetterlandAll <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0,]#& Hydrolakeid_candidate$TypeID==1
GLCPdata_filtered_wetterlandAll  <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_wetterlandAll$Hylak_id, ]
#Hydrolakeid_wetterlandbothimpacts <- Hydrolakeid_candidate[Hydrolakeid_candidate$Unuse==0 & Hydrolakeid_candidate$WSlmPvalue<0.05 & Hydrolakeid_candidate$WSlmRate<0 & Hydrolakeid_candidate$TypeID==0,]#& Hydrolakeid_candidate$TypeID==1
#GLCPdata_filtered_wetterlandbothimpacts <- GLCPdata_filtered[GLCPdata_filtered$Hylak_id %in% Hydrolakeid_wetterlandbothimpacts$Hylak_id, ]

sum(GLCPdata_filtered_wetterlandfreshwater$pop_sum)
sum(GLCPdata_filtered_wetterlandenvidegradation$pop_sum)
sum(GLCPdata_filtered_wetterlandenergyreduction$pop_sum)
sum(GLCPdata_filtered_humidlandAll$pop_sum)


#Figure plot
##########bg20210507
##########Figure 5 affected population living in lake basins with significant changes
specie <-c (rep("0Global",4),rep("1Dry land",4), rep("2Wet land",4))
condition <- rep(c("0 Freshwater decline","1 Environmental degradation","2 Either"),3)
#in lake basins

specie2 <-c (rep("0Global",4),rep("1Arid land",4), rep("2Humid land",4))
condition2 <- rep(c("0 Freshwater decline","1 Environmental degradation","2 Energy reduction","3 Any"),3)

value2 <- c(sum(GLCPdata_filtered_globalfreshwater$pop_sum),sum(GLCPdata_filtered_envidegradation$pop_sum),sum(GLCPdata_filtered_energyreduction$pop_sum),sum(GLCPdata_filtered_globalAll$pop_sum),
           sum(GLCPdata_filtered_drylandfreshwater$pop_sum),sum(GLCPdata_filtered_drylandenvidegradation$pop_sum),sum(GLCPdata_filtered_drylandenergyreduction$pop_sum),sum(GLCPdata_filtered_drylandAll$pop_sum),
           sum(GLCPdata_filtered_wetterlandfreshwater$pop_sum),sum(GLCPdata_filtered_wetterlandenvidegradation$pop_sum),sum(GLCPdata_filtered_wetterlandenergyreduction$pop_sum),sum(GLCPdata_filtered_wetterlandAll$pop_sum))/1000000000
data2 <- data.frame(specie2,condition2,value2)

#Record 210830
#Dry: 424334502, 312482458, 111844429
#Wet: 1031898227, 338704639, 328854258

#2553667480.927961, 
#7300000000
Gloal.pop=7300000000/100
Dry.pop=2536689804.333329/100
Wet.pop=Gloal.pop-Dry.pop

value3 <- c(sum(GLCPdata_filtered_globalfreshwater$pop_sum)/Gloal.pop,sum(GLCPdata_filtered_envidegradation$pop_sum)/Gloal.pop,sum(GLCPdata_filtered_energyreduction$pop_sum)/Gloal.pop,sum(GLCPdata_filtered_globalAll$pop_sum)/Gloal.pop,
            sum(GLCPdata_filtered_drylandfreshwater$pop_sum)/Dry.pop,sum(GLCPdata_filtered_drylandenvidegradation$pop_sum)/Dry.pop,sum(GLCPdata_filtered_drylandenergyreduction$pop_sum)/Dry.pop,sum(GLCPdata_filtered_drylandAll$pop_sum)/Dry.pop,
            sum(GLCPdata_filtered_wetterlandfreshwater$pop_sum)/Wet.pop,sum(GLCPdata_filtered_wetterlandenvidegradation$pop_sum)/Wet.pop,sum(GLCPdata_filtered_wetterlandenergyreduction$pop_sum)/Wet.pop,sum(GLCPdata_filtered_wetterlandAll$pop_sum)/Wet.pop)
data3 <- data.frame(specie2,condition2,value3)

#specie4 <- c(rep("0Global",4),rep("1Arid land",4), rep("2Humid land",4))  
#condition4 <- rep(c("0 Freshwater decline","1 Environmental degradation","2 Energy reduction","3 Any"),3)

specie4 <- c(rep("0Freshwater decline",3),rep("1 Environmental degradation",3), rep("2 Energy reduction",3), rep("3 Any",3))  
condition4 <- rep(c("0Global","1Arid regions","2 Humid regions"),4)
value4 <- c(sum(GLCPdata_filtered_globalfreshwater$pop_sum)/Gloal.pop,sum(GLCPdata_filtered_drylandfreshwater$pop_sum)/Dry.pop,sum(GLCPdata_filtered_wetterlandfreshwater$pop_sum)/Wet.pop,
            sum(GLCPdata_filtered_envidegradation$pop_sum)/Gloal.pop,sum(GLCPdata_filtered_drylandenvidegradation$pop_sum)/Dry.pop,sum(GLCPdata_filtered_wetterlandenvidegradation$pop_sum)/Wet.pop,
            sum(GLCPdata_filtered_energyreduction$pop_sum)/Gloal.pop,sum(GLCPdata_filtered_drylandenergyreduction$pop_sum)/Dry.pop,sum(GLCPdata_filtered_wetterlandenergyreduction$pop_sum)/Wet.pop,
            sum(GLCPdata_filtered_globalAll$pop_sum)/Gloal.pop,sum(GLCPdata_filtered_drylandAll$pop_sum)/Dry.pop,sum(GLCPdata_filtered_wetterlandAll$pop_sum)/Wet.pop)
data4 <- data.frame(specie4,condition4,value4)

#Same but in percentage, FinalUsed 211116
jpeg(Outfigure, 
     width = 6.4, height = 3.6, units = 'in', res = 600)
ggplot(data4, aes(x=as.factor(specie4), fill=as.factor(condition4),y=value4, width =0.75)) +  
  geom_bar(position="dodge",stat = "identity", show.legend = TRUE) +
  geom_hline(yintercept=0)+
  #geom_text(aes(x = as.factor(specie), y = yvalues + 0.9*sign(yvalues), label = round(yvalues, 2)))+
  scale_fill_manual(values = alpha(c("#dedad8","#c85856","#6774d1", "#6ec3c1","#0d5f8a","#9dcc5f","#335120","#6ec3c1","#0d5f8a","#9dcc5f","#335120","#6ec3c1"),1),name = "", labels=c("Global","Arid regions","Humid regions")) +
  #scale_fill_manual(labels=c("Freshwater decline","Environmental degradation","Either"))+
  #scale_fill_viridis(discrete = T)+
  labs(title="",y=bquote('Affected population (%)'), x="",size=12) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,29))+
  scale_y_continuous(expand = c(0,0),limits = c(0,29),breaks = seq(0,25, by=5))+
  scale_x_discrete(labels=c("Freshwater\ndecline","Environmental\ndegradation","Hydropower\nenergy reduction","Any"))+
  #scale_fill_discrete(name = "", labels=c("Water loss","Water gain"))+
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends="last", type = "closed")),
                     axis.text.x = element_text(angle = 0, vjust = 1, hjust=1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) 
dev.off()




global.fresh <- (sum(GLCPdata_filtered_globalfreshwater$pop_sum)-sum(GLCPdata_filtered_bothimpacts$pop_sum))/Gloal.pop
global.envi <- (sum(GLCPdata_filtered_envidegradation$pop_sum)-sum(GLCPdata_filtered_bothimpacts$pop_sum))/Gloal.pop
global.both <- sum(GLCPdata_filtered_bothimpacts$pop_sum)/Gloal.pop
global.all <- global.fresh+global.envi+global.both
dry.fresh <- (sum(GLCPdata_filtered_drylandfreshwater$pop_sum)-sum(GLCPdata_filtered_drylandbothimpacts$pop_sum))/Dry.pop
dry.envi <- (sum(GLCPdata_filtered_drylandenvidegradation$pop_sum)-sum(GLCPdata_filtered_drylandbothimpacts$pop_sum))/Dry.pop
dry.both <- sum(GLCPdata_filtered_drylandbothimpacts$pop_sum)/Dry.pop
dry.all <- dry.fresh+dry.envi+dry.both
wet.fresh <- (sum(GLCPdata_filtered_wetterlandfreshwater$pop_sum)-sum(GLCPdata_filtered_wetterlandbothimpacts$pop_sum))/Wet.pop
wet.envi <- (sum(GLCPdata_filtered_wetterlandenvidegradation$pop_sum)-sum(GLCPdata_filtered_wetterlandbothimpacts$pop_sum))/Wet.pop
wet.both <- sum(GLCPdata_filtered_wetterlandbothimpacts$pop_sum)/Wet.pop
wet.all <- wet.fresh+wet.envi+wet.both

specie.separate <-c (rep("0Global",4),rep("1Dry land",4), rep("2Wet land",4))
condition.separate <- rep(c("0 Freshwater decline","1 Environmental degradation","2 Both","3 All"),3)
value.separate <- c(global.fresh,global.envi,global.both,global.all,
                    dry.fresh,dry.envi,dry.both,dry.all,
                    wet.fresh,wet.envi,wet.both,wet.all)
data.separate <- data.frame(specie.separate,condition.separate,value.separate)

#U210830 using percentage to highlight dryland impacts

#Grouped
#ggplot(data, aes(fill=condition, y=value, x=specie)) + geom_bar(position="dodge",stat = "identity")
#colour=c("blue","red", "blue", "red"),
#####Either
jpeg(Outfigure, 
     width = 6.4, height = 3.6, units = 'in', res = 600)
ggplot(data2, aes(x=as.factor(specie2), fill=as.factor(condition2),y=value2, width =0.75)) +  
  geom_bar(position="dodge",stat = "identity", show.legend = TRUE) +
  geom_hline(yintercept=0)+
  #geom_text(aes(x = as.factor(specie), y = yvalues + 0.9*sign(yvalues), label = round(yvalues, 2)))+
  scale_fill_manual(values = alpha(c("#0d5f8a","#9dcc5f","#335120", "#6ec3c1","#0d5f8a","#9dcc5f","#335120","#6ec3c1","#0d5f8a","#9dcc5f","#335120","#6ec3c1"),1),name = "", labels=c("Freshwater decline","Environmental degradation","Energy reduction (hydropower)","Any of the above")) +
  #scale_fill_manual(labels=c("Freshwater decline","Environmental degradation","Either"))+
  #scale_fill_viridis(discrete = T)+
  labs(title="",y=bquote('Number of affected people \n (billions)'), x="",size=12) +
  scale_y_continuous(expand = c(0,0),limits = c(0,2.3))+
  scale_x_discrete(labels=c("Global","Arid land","Humid land"))+
  #scale_fill_discrete(name = "", labels=c("Water loss","Water gain"))+
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends="last", type = "closed")),
                     axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) 
dev.off()
#Same but in percentage
jpeg(Outfigure, 
     width = 6.4, height = 3.6, units = 'in', res = 600)
ggplot(data3, aes(x=as.factor(specie2), fill=as.factor(condition2),y=value3, width =0.75)) +  
  geom_bar(position="dodge",stat = "identity", show.legend = TRUE) +
  geom_hline(yintercept=0)+
  #geom_text(aes(x = as.factor(specie), y = yvalues + 0.9*sign(yvalues), label = round(yvalues, 2)))+
  scale_fill_manual(values = alpha(c("#0d5f8a","#9dcc5f","#335120", "#6ec3c1","#0d5f8a","#9dcc5f","#335120","#6ec3c1","#0d5f8a","#9dcc5f","#335120","#6ec3c1"),1),name = "", labels=c("Freshwater decline","Environmental degradation","Hydropower energy reduction","Any of the above")) +
  #scale_fill_manual(labels=c("Freshwater decline","Environmental degradation","Either"))+
  #scale_fill_viridis(discrete = T)+
  labs(title="",y=bquote('Affected global population (%)'), x="",size=12) +
  #scale_y_continuous(expand = c(0,0),limits = c(0,29))+
  scale_y_continuous(expand = c(0,0),limits = c(0,29),breaks = seq(0,25, by=5))+
  scale_x_discrete(labels=c("Global","Arid regions","Humid regions"))+
  #scale_fill_discrete(name = "", labels=c("Water loss","Water gain"))+
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends="last", type = "closed")),
                     axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) 
dev.off()


######Both
# jpeg(Outfigure, 
#      width = 5.4, height = 3.6, units = 'in', res = 600)
# ggplot(data2, aes(x=as.factor(specie2), fill=as.factor(condition2),y=value2, width =0.75)) +  
#   geom_bar(position="dodge",stat = "identity", show.legend = TRUE) +
#   geom_hline(yintercept=0)+
#   #geom_text(aes(x = as.factor(specie), y = yvalues + 0.9*sign(yvalues), label = round(yvalues, 2)))+
#   scale_fill_manual(values = alpha(c("#340042","#1f8079","#fce41e", "#340042","#1f8079","#fce41e","#340042","#1f8079","#fce41e"),1),name = "", labels=c("Freshwater decline","Environmental degradation","Both")) +
#   #scale_fill_manual(labels=c("Freshwater decline","Environmental degradation","Either"))+
#   #scale_fill_viridis(discrete = T)+
#   labs(title="",y=bquote('Number of affected people \n (billions)'), x="",size=12) +
#   scale_y_continuous(expand = c(0,0),limits = c(0,2))+
#   scale_x_discrete(labels=c("Global","Dry land","Wet land"))+
#   #scale_fill_discrete(name = "", labels=c("Water loss","Water gain"))+
#   theme_bw() + theme(axis.line = element_line(colour = "black"),
#                      axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends="last", type = "closed")),
#                      axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
#                      panel.grid.major = element_blank(),
#                      panel.grid.minor = element_blank(),
#                      panel.border = element_blank(),
#                      panel.background = element_blank()) 
# dev.off()




jpeg(Outfigure, 
     width = 5.4, height = 3.6, units = 'in', res = 600)
ggplot(data.separate, aes(x=as.factor(specie.separate), fill=as.factor(condition.separate),y=value.separate, width =0.75)) +  
  geom_bar(position="dodge",stat = "identity", show.legend = TRUE) +
  geom_hline(yintercept=0)+
  #geom_text(aes(x = as.factor(specie), y = yvalues + 0.9*sign(yvalues), label = round(yvalues, 2)))+
  scale_fill_manual(values = alpha(c("#340042","#1f8079","#fce41e","3e3e3e", "#340042","#1f8079","#fce41e","3e3e3e","#340042","#1f8079","#fce41e","3e3e3e"),1),name = "", labels=c("Freshwater decline","Environmental degradation","Both","Sum of above")) +
  #scale_fill_manual(labels=c("Freshwater decline","Environmental degradation","Either"))+
  #scale_fill_viridis(discrete = T)+
  labs(title="",y=bquote('Percentage of affected people (%)'), x="",size=12) +
  scale_y_continuous(expand = c(0,0),limits = c(0,25))+
  scale_x_discrete(labels=c("Global","Dry land","Wet land"))+
  #scale_fill_discrete(name = "", labels=c("Water loss","Water gain"))+
  theme_bw() + theme(axis.line = element_line(colour = "black"),
                     axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"),ends="last", type = "closed")),
                     axis.text.x = element_text(angle = 40, vjust = 1, hjust=1),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()) 
dev.off()

TendR2csv <- '/Users/fangfangyao/Documents/Proj/GlakesVolume/GLCP/Input/Glakes_SignificantTrends.csv'
Outfigure2 < "/Users/fangfangyao/Documents/Proj/GlakesVolume/GLCP/Output/Figure5b_SignificantTrends.jpg"
TendR2data <- read.csv(file=TendR2csv)

ggplot(TendR2data, aes(x=InArid, y=TrendR2)) + 
  geom_boxplot()+
  labs(title="",y=bquote('R[2] values of the regressed linear trends'), x="",size=12) +
  theme_bw() 


p <-ggplot(mydata,aes(x_values, values))
p +geom_bar(stat = "identity", aes(fill = type), position = "dodge") +
  xlab("Monthly Coverage") + ylab("Count of Reservoirs") +
  ggtitle("") +
  theme_bw()
ggsave(outputfilefull,width = 8, height = 6,units='in', dpi=300)
axis.title.y = element_text(size = rel(1.5), angle = 90)
