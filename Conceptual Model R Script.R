#Amanda C. Castaing 
#Final Project 
#Transport of root produced Glucose and Microbial Biomass diffusion through the soil starting at the root surface. 
#Soil Microbial Dynamics & Glucose Diffusion Model 

  #Load in libraries ----
library(tidyverse); library(RColorBrewer); library(deSolve)

  # Function: Diffusion for Glucose and Microbes
  RImodel <- function(time, state, parms, distance, nr0) {
    with(as.list(parms), {
      nr = nr0 #renaming this argument, throws error of 'nr' is used as argument name in ode.1D
      dr = distance/nr
      r <- seq(dr/2,distance - dr/2,dr) #all the midpoints/observations
      
      GLUCOSE <- state[1:nr]
      MICROBE <- state[(nr+1):(2*nr)]
      
      #Glucose Diffusion Equation
      d1 = GLUCOSE[1:(nr-2)] - 2*GLUCOSE[2:(nr-1)] + GLUCOSE[3:nr]
      d2 = GLUCOSE[3:nr] - GLUCOSE[1:(nr-2)]
      dGLUCOSE = (d1/dr**2 + d2/2/dr/r[2:(nr-1)]) * D
      
      
      #Boundary Conditions
      b1 <- D/(dr^2) * (-GLUCOSE[1] + GLUCOSE[2]) + D/(2*r[1]*dr) * (GLUCOSE[2] - GLUCOSE[1])
      b2 <- D/(dr^2) * (-GLUCOSE[nr] + GLUCOSE[nr-1]) + D/(2*r[nr]*dr) * (GLUCOSE[nr] - GLUCOSE[nr-1])
      
      #Derivative vector 
      dGLUCOSE = c(b1, dGLUCOSE, b2)
      
  # Biology: Glucose and Microbes
      GlucoseUptake <- uptakeRate * GLUCOSE * MICROBE 
      MicrobeGrowth <- GlucoseUptake * CUE 
      MicrobeMortality <- mortalityRate * MICROBE
      
  # Rate of change = Flux gradient + Biology
      dGLUCOSE <-  dGLUCOSE - GlucoseUptake + MicrobeMortality
      dMICROBE <- MicrobeGrowth - MicrobeMortality
      
      return(list(c(dGLUCOSE, dMICROBE)))
    })
  }

#Time
dt = 3600/4 #time step (seconds)
days = 10 #number of days to run the simulation
duration = days*86400 #days to seconds
nt = duration/dt +1 #number of observations
time = seq(0,duration,dt) #seconds

#Space (1-D) 
r0 = 0.03 #where root meets soil (mm)
distance = 30 #distance away from center (mm) 

#Parameters
area_volume <-  pi * distance^2 #mm^3
#Since we are working in 1-D, area = volume of a sub-cylinder 
#sourced from Leady et al. 1997
bulk_density <- 1.01 #g/cm^3 = mg/mm^3
#Gainesville average according to Hagan et al (2021)
soil_mass <- bulk_density * area_volume #mg
SOC <- soil_mass * 0.02 #Assumes SOC is 2% (mg)

C_add <- SOC * 0.01 #Add 1% of SOC in terms of glucose carbon (mg)
glu_needed <- C_add * 180.156 / 72.06 * 1000 #glucose plant production (mg)

D <- 0.5 * 10^-4 #diffusion of glucose mm2/s #sourced from Chenu and Roberson (1996) given VWC ~ 30

# Model application

# Model parameters
R <- 20
N <- 100
dr <- R/N
r <- seq(dr/2, by = dr, len = N)
ri <- seq(0, by = dr, len = N+1)
dri <- dr
parms <- c(D = 0.5 * 10^-4,
           uptakeRate = 0.2/86400,
           mortalityRate = 0.2/86400,
           CUE = 0.5)


# Initial conditions
Glucose = rep(0, N)
Glucose[1] <- C_add <- SOC * 0.01 
Microbes = rep(C_add/100, N)

state = c(Glucose,Microbes)

# Time points for output


# Solve the model using the default banded reformulation method
model <- RImodel(0, state, parms, distance, N)
result <- ode.1D(y = state, times = time, func = RImodel, parms = parms, distance = R, nr0 = N, dimens = N,
               nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time = result[,1]

# Extracting results
GLUCOSE_result <- result[, 2:(N + 1)]
MICROBE_result <- result[, (N + 2):(2*N + 1)]


#----------------------------------------------
#Plots Possibilities

#Plot Glucose concentration and Microbe population over time
plot(time, GLUCOSE_result[, 1], type = "l", col = "orange", xlab = "Time", ylab = "Glucose Concentration", main = "Glucose Concentration and Microbe Population Over Time at [,1]")
lines(time, MICROBE_result[, 1], type = "l", col = "green")
legend("topright", legend = c("Glucose", "Microbe"), col = c("orange", "green"), lty = 1)

plot(time, GLUCOSE_result[, 50], type = "l", col = "orange", xlab = "Time", ylab = "Glucose Concentration", main = "Glucose Concentration and Microbe Population Over Time at [,50]")
lines(time, MICROBE_result[, 50], type = "l", col = "green")
legend("topright", legend = c("Glucose", "Microbe"), col = c("orange", "green"), lty = 1)


plot(time, GLUCOSE_result[, 100], type = "l", col = "orange", xlab = "Time", ylab = "Glucose Concentration", main = "Glucose Concentration and Microbe Population Over Time at [,100]")
lines(time, MICROBE_result[, 100], type = "l", col = "green")
legend("topright", legend = c("Glucose", "Microbe"), col = c("orange", "green"), lty = 1)




#Calculate the percentage of glucose and microbes at different locations
percentage_glucose <- (GLUCOSE_result / mean(GLUCOSE_result)) * 100
percentage_microbes <- (MICROBE_result / mean(MICROBE_result)) * 100


# Plot the percentages
par(mfrow = c(2, 1)) # Create a 2x1 grid of plots

# Percentage of Glucose at Different Locations
plot(r, percentage_glucose[1,], type = "l", 
     xlab = "Distance from Root (mm)",
     ylab = "Percentage",
     main = paste("Percentage of Glucose\nTime:", time[1]))

for (i in 1:4) {
  lines(r, percentage_glucose[i, ], type = "l",
        col= i)
      
}

plot(r, log10(percentage_glucose[1,]), type = "l", 
     xlab = "Distance from Root (mm)",
     ylab = "Percentage",
     main = paste("Percentage of Glucose\nTime:", time[1]))

for (i in 1:4) {
  lines(r, log10(percentage_glucose[i, ]), type = "l",
        col= i)
 
}


# Percentage of Microbe at different locations

# Plot the percentages
par(mfrow = c(2, 1)) # Create a 2x1 grid of plots


plot(r, percentage_microbes[1,], type = "l", 
     xlab = "Distance from Root (mm)",
     ylab = "Percentage",
     main = paste("Percentage of Microbe\nTime:", time[1]))


for (i in 1:4) {
  lines(r, percentage_microbes[i, ], type = "l",
        col= i)
 
}

plot(r, log10(percentage_microbes[1,]), type = "l", 
     xlab = "Distance from Root (mm)",
     ylab = "Percentage",
     main = paste("Percentage of Microbes\nTime:", time[1]))

for (i in 1:4) {
  lines(r, log10(percentage_microbes[i, ]), type = "l",
        col= i)
  
}


#Calculate the distance between glucose and microbes at different locations
distance_between <- abs(percentage_glucose - percentage_microbes)

# Plot the distance between glucose and microbes
par(mfrow = c(2, 1)) # Create a 2x1 grid of plots

# Distance between Glucose and Microbes at Different Locations
for (i in 1:4) {
  plot(r, distance_between[i,], type = "l", 
       col = "brown",
       xlab = "Distance from Root (mm)", 
       ylab = "Distance",
       main = paste("Distance Between Glucose and Microbes\nTime:", time[i]))
}

par(mfrow = c(1, 1)) # Reset plotting parameters



#Glucose Diffusion Profile Over Time
par(mfrow = c(2, 2)) # Create a 2x2 grid of plots
for (i in 1:4){
  plot(r, GLUCOSE_result[i*25,] , type = "l", 
       col = "orange",
       xlab = "Distance from Root (mm)", 
       ylab = "Glucose Concentration",
       main = paste("Glucose Diffusion Profile\nTime:", time[i*25]))
}

#Microbial Population Growth Over Time
plot(time, MICROBE_result[, N], type = "l", 
     col = "green",
     xlab = "Time", 
     ylab = "Microbial Population",
     main = "Microbial Population Growth Over Time")

#Glucose Uptake by Microbes
GlucoseUptake <- parms["uptakeRate"] * GLUCOSE_result * MICROBE_result
plot(time, GlucoseUptake[, N], type = "l", 
     col = "yellow",
     xlab = "Time", 
     ylab = "Glucose Uptake Rate",
     main = "Glucose Uptake by Microbes")

#Microbial Growth and Mortality Rates
MicrobeGrowth <- parms["uptakeRate"] * GLUCOSE_result * MICROBE_result * parms["CUE"]
MicrobeMortality <- parms["mortalityRate"] * MICROBE_result

plot(time, MicrobeGrowth[, N], type = "l", 
     xlab = "Time", ylab = "Rate",
     main = "Microbial Growth and Mortality Rates")
lines(time, MicrobeMortality[, N], col = "red")
legend("topright", legend = c("Growth", "Mortality"), col = c("black", "red"), lty = 1)

#Glucose Concentration and Microbial Population Correlation
plot(GLUCOSE_result[length(time),], MICROBE_result[length(time),], col = "darkturquoise", pch = 16,
     xlab = "Glucose Concentration", 
     ylab = "Microbial Population",
     main = "Glucose-Microbial Correlation")

#Total Microbial Biomass Over Time
total_microbial_biomass <- rowSums(MICROBE_result)
plot(time, total_microbial_biomass, type = "l", 
     col = "purple",
     xlab = "Time",
     ylab = "Total Microbial Biomass",
     main = "Total Microbial Biomass Over Time")

#-------------------------------
# Comapring Parameters higher & low changes in Value 

parms <- c(D = 0.75 * 10^-4,
          uptakeRate = 0.3/86400,
          mortalityRate = 0.3/86400,
          CUE = 0.7)

par(mfrow = c(1, 1)) # Reset plotting parameters



#Glucose Diffusion Profile Over Time
par(mfrow = c(2, 2)) # Create a 2x2 grid of plots
for (i in 1:4){
  plot(r, GLUCOSE_result[i*25,] , type = "l", 
       col = "orange",
       xlab = "Distance from Root (mm)", 
       ylab = "Glucose Concentration",
       main = paste("Glucose Diffusion Profile\nTime:", time[i*25]))
}

#Microbial Population Growth Over Time
plot(time, MICROBE_result[, N], type = "l", 
     col = "green",
     xlab = "Time", 
     ylab = "Microbial Population",
     main = "Microbial Population Growth Over Time")

#Glucose Uptake by Microbes
GlucoseUptake <- parms["uptakeRate"] * GLUCOSE_result * MICROBE_result
plot(time, GlucoseUptake[, N], type = "l", 
     col = "yellow",
     xlab = "Time", 
     ylab = "Glucose Uptake Rate",
     main = "Glucose Uptake by Microbes")

#Microbial Growth and Mortality Rates
MicrobeGrowth <- parms["uptakeRate"] * GLUCOSE_result * MICROBE_result * parms["CUE"]
MicrobeMortality <- parms["mortalityRate"] * MICROBE_result

plot(time, MicrobeGrowth[, N], type = "l", 
     xlab = "Time", ylab = "Rate",
     main = "Microbial Growth and Mortality Rates")
lines(time, MicrobeMortality[, N], col = "red")
legend("topright", legend = c("Growth", "Mortality"), col = c("black", "red"), lty = 1)

#Glucose Concentration and Microbial Population Correlation
plot(GLUCOSE_result[length(time),], MICROBE_result[length(time),], col = "darkturquoise", pch = 16,
     xlab = "Glucose Concentration", 
     ylab = "Microbial Population",
     main = "Glucose-Microbial Correlation")

#Total Microbial Biomass Over Time
total_microbial_biomass <- rowSums(MICROBE_result)
plot(time, total_microbial_biomass, type = "l", 
     col = "purple",
     xlab = "Time",
     ylab = "Total Microbial Biomass",
     main = "Total Microbial Biomass Over Time")


#--------------------------------------------------
parms <- c(D = 0.25 * 10^-4,
           uptakeRate = 0.1/86400,
           mortalityRate = 0.1/86400,
           CUE = 0.3)
par(mfrow = c(1, 1)) # Reset plotting parameters



#Glucose Diffusion Profile Over Time
par(mfrow = c(2, 2)) # Create a 2x2 grid of plots
for (i in 1:4){
  plot(r, GLUCOSE_result[i*25,] , type = "l", 
       col = "orange",
       xlab = "Distance from Root (mm)", 
       ylab = "Glucose Concentration",
       main = paste("Glucose Diffusion Profile\nTime:", time[i*25]))
}

#Microbial Population Growth Over Time
plot(time, MICROBE_result[, N], type = "l", 
     col = "green",
     xlab = "Time", 
     ylab = "Microbial Population",
     main = "Microbial Population Growth Over Time")

#Glucose Uptake by Microbes
GlucoseUptake <- parms["uptakeRate"] * GLUCOSE_result * MICROBE_result
plot(time, GlucoseUptake[, N], type = "l", 
     col = "yellow",
     xlab = "Time", 
     ylab = "Glucose Uptake Rate",
     main = "Glucose Uptake by Microbes")

#Microbial Growth and Mortality Rates
MicrobeGrowth <- parms["uptakeRate"] * GLUCOSE_result * MICROBE_result * parms["CUE"]
MicrobeMortality <- parms["mortalityRate"] * MICROBE_result

plot(time, MicrobeGrowth[, N], type = "l", 
     xlab = "Time", ylab = "Rate",
     main = "Microbial Growth and Mortality Rates")
lines(time, MicrobeMortality[, N], col = "red")
legend("topright", legend = c("Growth", "Mortality"), col = c("black", "red"), lty = 1)

#Glucose Concentration and Microbial Population Correlation
plot(GLUCOSE_result[length(time),], MICROBE_result[length(time),], col = "darkturquoise", pch = 16,
     xlab = "Glucose Concentration", 
     ylab = "Microbial Population",
     main = "Glucose-Microbial Correlation")

#Total Microbial Biomass Over Time
total_microbial_biomass <- rowSums(MICROBE_result)
plot(time, total_microbial_biomass, type = "l", 
     col = "purple",
     xlab = "Time",
     ylab = "Total Microbial Biomass",
     main = "Total Microbial Biomass Over Time")
