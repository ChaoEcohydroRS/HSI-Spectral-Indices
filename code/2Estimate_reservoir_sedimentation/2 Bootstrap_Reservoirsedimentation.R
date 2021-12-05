###################################################################################################################
##Bootstrapping for estimating the mean sedimentation rate of large reservoirs and the confidence interval#########
##bootstrapping refer to Efron, B. & Tibshirani, R. J. An introduction to the bootstrap. (CRC press, 1994)#########
###################################################################################################################
setwd("G:/My Drive/0Desktop/Program/R/GLakes30m/GlakesVolumeCleaned/Glakes/")

library(MASS)
library(sm)
library(ggplot2)

##input: number of iterations for bootstrapping 
num.iteration = 1500
##input: in situ sedimentation survey data from United States Army Corps of Engineers
Sedimnetsurvey = "./data/Preprocessed/Sedimentation_Survery/USArmyCorps_reservoir_sedimentation_annualrate.txt"
##output: distribution of the simulate mean sedimentation rate
Outjpeg = "./data/Output_sedimentation/Reservoir_Sedimentation_meanrate_Bootstrap.jpg"

data=matrix(scan(Sedimnetsurvey),ncol=2,byrow=T)

xdata=data[,2]


#suppose the data is in the vector 'xdata'
N=length(xdata)		#number of data points.

nsim=num.iteration	#number of simulations, each of length N

xeval=seq(min(xdata)-sd(xdata),max(xdata)+sd(xdata),length=100)
neval=length(xeval)
xsimpdf=matrix(0,nrow=nsim,ncol=neval)

#initialize vectors for the statistics
meansim=1:nsim
mediansim=1:nsim
sdsim=1:nsim
skewsim=1:nsim
iqrsim=1:nsim

# Evaluate the Kernel density PDF at the evaluation points based on the data
xdensityorig=sm.density(xdata, eval.points=xeval, display="none")$estimate

#compute bandwidth of the data
band=hnorm(xdata)
######## Simulation ##########


for(i in 1:nsim){
    
    # Simulate from the desired model (i.e., PDF) 
    # each simulation should of the same length as the original data
    # used to fit the model, N = length(xdata)
    
    #simulate from empirical CDF - BOOTSTRAP
    N = length(xdata)
    xsim=sample(xdata,N, replace=TRUE)
    
    #From Kernel PDF  - smooth bootstrap..
    #xsim=sample(xdata,N, replace=TRUE) + band*rnorm(N)
    
    #From Kernel PDF - Smooth bootstrap with correction to reproduce the variance
    
    #varkern is the variance of the kernel you choose..
    #varkern=1		#variance of Normal kernel is 1
    #denom=sqrt(1 + (band*band*varkern/var(xdata)))
    #xsim=mean(xdata) + ((sample(xdata,N,replace=TRUE) - mean(xdata) + band*rnorm(N))/denom)
    
    # compute the statistics from the simulation
    meansim[i]=mean(xsim)
    sdsim[i]=sd(xsim)
    #sdsim[i]=sd(log(xsim))
    #skewsim[i]=skew(xsim)
    iqrsim[i]=diff(quantile(xsim,c(0.25,0.75)))
    mediansim[i]=quantile(xsim,c(0.5))
    
    #estimate the Kernel PDF at the evaluation points based on the simulated data
    
    xsimpdf[i,]=sm.density(xsim, eval.points=xeval, display="none")$estimate
    
}

xlabel <-c (rep("Reservoirs",nsim))
data<-data.frame(Dryland=xlabel,SedimentAR=meansim)
jpeg(Outjpeg, 
     width = 3.5, height = 4, units = 'in', res = 600)
ggplot(data, aes(x=Dryland,y=SedimentAR,color=Dryland)) +
    geom_boxplot(width=0.4,fill='#A4A4A4', color="black")+
    labs(title="",y=bquote('Mean rate of annual storage loss (%)'), x="",size=16)+
    theme_classic()

dev.off()
#quantile(meansim, c(.02, .50, .97)) 

