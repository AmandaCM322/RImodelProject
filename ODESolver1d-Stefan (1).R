#Shelby Kucharski
#University of Florida | Agronomy
#August 2022 (last revised 09/06/2022)

#Transport of root produced glucose/carbon through soil - Final project using ODE

#Load in libraries ----
library(tidyverse); library(RColorBrewer); library(deSolve)

#Function ----
RImodel <- function(time, conc, D, distance, nr0, dr){
  nr = nr0 #renaming this argument, throws error of 'nr' is used as argument name in ode.1D
  r <- seq(dr/2,distance - dr/2,dr) #all the midpoints/observations
  
  d1 <- conc[1:(nr-2)] - 2*conc[2:(nr-1)] + conc[3:nr]
  d2 <- conc[3:nr] -conc[1:(nr-2)]
  dconc <- (d1/dr^2 + d2/2/dr/r[2:(nr-1)])*D 
  
  #boundary conditions
  b1 = D/(dr^2) * (-conc[1] + conc[2]) + D/(2*r[1]*dr) *  (conc[2]-conc[1])
  b2 = D/(dr^2) * (-conc[nr]+ conc[nr-1]) + D/(2*r[nr]*dr) * (conc[nr]-conc[nr-1])
  
  #complete the derivative vector
  dconc <- c(b1,dconc,b2)
  
  
#  for (i in 2:(nr-1)){
#    conc[i] <- D/(dr^2) * (-2*conc[i] + conc[i-1] + conc[i+2]) + D/(2*r[i]*dr) * (conc[i+2] - conc[i-1])
#  }
#  conc[1] <- D/(dr^2) * (-conc[1] + conc[2]) + D/(2*r[1]*dr) *  (conc[2]-conc[1]) #no flux in 
#  conc[nr] <- D/(dr^2) * (-conc[nr]+ conc[nr-1]) + D/(2*r[nr]*dr) * (conc[nr]-conc[nr-1]) #flux out=0 at boundary
  #conc[,i+1] <- conc[,i] + dconc*dt 
  return(list(c(dconc)))
}

#Dimensions ----
#Time
dt <- 3600/4 #time step (seconds)
days <- 10 #number of days to run the simulation
duration <- days*86400 #days to seconds
nt <- duration/dt +1 #number of observations
time <- seq(0,duration,dt) #seconds

#Space (1-D) 
r0 <- 0 #where root meets soil (mm)
distance <- 30 #distance away from center (mm) 
nr <- distance #number of observations
dr <- distance/nr #distance between each observation (mm)

#Parameters ----
D <- 0.5 * 10^-4 #diffusion of glucose mm2/s 
  #sourced from Chenu and Roberson (1996) given VWC ~ 30

#Initial Condition ----
area_volume <-  pi * distance^2 #mm^3
  #Since we are working in 1-D, area = volume of a sub-cylinder 
  #sourced from Leady et al. 1997
bulk_density <- 1.01 #g/cm^3 = mg/mm^3
 #Gainesville average according to Hagan et al (2021)
soil_mass <- bulk_density * area_volume #mg
SOC <- soil_mass * 0.02 #Assumes SOC is 2% (mg)
C_add <- SOC * 0.01 #Add 1% of SOC in terms of glucose carbon (mg)
plant_production <- C_add #Carbon released via glucose addition (mg)
  #glucose plant production would be glu_needed <- C_add * 180.156 / 72.06 * 1000 (mg)
conc <- numeric(nr)
conc[1] <- plant_production/dr #concentration (mg/mm^3)

#ODE ----
model <- RImodel(time, conc, D, distance, nr, dr)
result <- ode.1D(y=conc, times=time, func=RImodel, parms=D, dimens=30, method="lsoda",distance=distance, nr0=nr,dr=dr)

r <- seq(dr/2,distance - dr/2,dr)
time=result[,1]
y = result[,2:(length(r)+1)]


plot(r,y[length(time),])


