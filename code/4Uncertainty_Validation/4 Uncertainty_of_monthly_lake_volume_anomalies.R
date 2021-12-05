###################################################################################################################
##Estimating errors of lake water storage via propagating errors in levels, areas, and hypsometric curve fitting###
##Refer to http://ipl.physics.harvard.edu/wp-uploads/2013/03/PS3_Error_Propagation_sp13.pdf for error propagation##
###################################################################################################################
setwd("G:/My Drive/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Glakes/")

library(XLConnect)
library(mblm)
library(ggplot2)
options("scipen"=100, "digits"=16)

##input: time series of lake areas, levels, and volumes in one csv file
inputcsv="./data/Output_Volume_TimeSeries/Monthly/ID82Large_aral_sea_West30m_GEE1992_2020_monthlyVolume_Poly.csv"
##output: lake water storage errors
outputcsv="./data/Output_Volume_TimeSeries/Uncertainty/ID82Large_aral_sea_West30m_GEE1992_2020_monthlyVolume_Error.csv"

##Get parameters of the fitted hypsometric curve
##Area=-78860948454.61244+8695784894.057913*x-314080884.7731609*x^2+3893584.76621067*x^3
kEval1=8695784894.057913
kEval2=-314080884.7731609
kEval3=3893584.76621067
Intercept=-78860948454.61244
H.RMSE=7.50642e+07   
Mean.Levelerr = 0.00920644 


##Mean error in mapped lake areas: 2.2%, obtained from Yao et al. 2019, Remote Sensing of Environment
err.mappedarea = 0.022
##Scale km2 to m2
Scalekm2_m2=1e6

mydata <- read.csv(file = inputcsv, header = TRUE, sep = ",")

## Get all available Level and area data
Level <- mydata$Mean_Levels
Level.err <- mydata$Mean_Level_Errors
Area <- mydata$Cleaned_Area*Scalekm2_m2
Time <- mydata$Time

mydata$Level.errOri <- mydata$Mean_Level_Errors

##Possible maximum err in level; customize if needed, larger number, longer computation time
levelerr.max = 10

##Combine errors in areas, levels (as provided by altimetry database), and hypsometry into level errors
for (index in 1:length(Level)) {   
    if(is.na(Level[index])) {
        next
    }
    area_currTimestep = 0
    if(is.na(Area[index]))
    {
        area_currTimestep=kEval1*Level[index]+kEval2*(Level[index]^2)+kEval3*(Level[index]^3)+Intercept
    }
    else
    {
        area_currTimestep=Area[index]
    }
    areaerr.comb <- sqrt(H.RMSE^2 + (area_currTimestep*err.mappedarea)^2)
    levelerr_currTimestep = 0 
    if(is.na(Level.err[index]))
    {
        levelerr_currTimestep=Mean.Levelerr
    }
    else
    {
        levelerr_currTimestep=Level.err[index]
    }
    conditionvalue.4loop = 1e9
    if(conditionvalue.4loop<areaerr.comb)
        conditionvalue.4loop=areaerr.comb+0.1
    for(levelerrtemp in seq(from=0, to=levelerr.max, by=0.001))
    {
        equation.x1 <-  kEval1*Level[index] - abs(min(Level,na.rm = TRUE)*kEval1)
        equation.x1.err <- abs(levelerrtemp*equation.x1/Level[index])
        equation.x2 <- kEval2*Level[index]^2 - abs(min(Level^2,na.rm = TRUE)*kEval2)
        equation.x2.err <- abs(2*levelerrtemp*equation.x2/Level[index])
        equation.x3 <-  kEval3*Level[index]^3 - abs(min(Level^3,na.rm = TRUE)*kEval3)
        equation.x3.err <- abs(3*levelerrtemp*equation.x3/Level[index])
        areaerr.simu <- sqrt(equation.x1.err^2+equation.x2.err^2+equation.x3.err^2)
        if(abs(areaerr.simu-areaerr.comb)<conditionvalue.4loop){
            
            mydata$Mean_Level_Errors[index] = sqrt(levelerrtemp^2+levelerr_currTimestep^2)
            Level.err[index]=mydata$Mean_Level_Errors[index]
            conditionvalue.4loop = abs(areaerr.simu-areaerr.comb)
        }
    }
}
mydata$Volume <- mydata$rws
meanVolume <- mean(mydata$Volume[!is.na(mydata$Volume)])
mydata$Volume[!is.na(mydata$Volume)]=mydata$Volume[!is.na(mydata$Volume)]-meanVolume

mydata$Volumeerr <- rep(NA,length(mydata$Volume))
mydata$TimeNum <- rep(1992.75,length(mydata$Volume))
mydata$Levelerr <- Level.err
mydata$Area <- mydata$Cleaned_Area

mydata$Volumeerrmax <- rep(NA,length(mydata$Volume))
##Estimate errors of lake water storage
for (index in 1:length(Level)) {  
    mydata$TimeNum[index] = mydata$TimeNum +index/12
    if(!is.na(Level[index]))
    {
        equation.x0 <- Intercept*Level[index] - min(Level,na.rm = TRUE)*Intercept
        equation.x0.err <- abs(Level.err[index]*equation.x0/Level[index])
        equation.x1 <-  kEval1*(Level[index]^2)/2 - min(Level^2,na.rm = TRUE)*kEval1/2
        equation.x1.err <- abs(2*Level.err[index]*equation.x1/Level[index])
        equation.x2 <- kEval2*(Level[index]^3)/3 - min(Level^3,na.rm = TRUE)*kEval2/3
        equation.x2.err <- abs(3*Level.err[index]*equation.x2/Level[index])
        equation.x3 <-  kEval3*(Level[index]^4)/4 - min(Level^4,na.rm = TRUE)*kEval3/4
        equation.x3.err <- abs(4*Level.err[index]*equation.x3/Level[index])
        volumeerr.comb <- sqrt(equation.x0.err^2+equation.x1.err^2+equation.x2.err^2+equation.x3.err^2)
        mydata$Volumeerr[index] = volumeerr.comb/1e9
    }
}
mydata[c("Time","TimeNum","Area","Volume","Volumeerr","Volumeerrmax","Levelerr")]
write.csv(mydata[c("Time","TimeNum","Area","Volume","Volumeerr","Volumeerrmax","Levelerr")],outputcsv)
    
ggplot(mydata, aes(x=TimeNum, y=Volume)) + 
    ggtitle("Estimated lake volume and volume errors") +
    xlab("Year") +
    ylab("Lake volume anomaly (Gt)") +
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=Volume-Volumeerr, ymax=Volume+Volumeerr), width=.2,position=position_dodge(0.05))


