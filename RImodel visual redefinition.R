#Load in libraries ----
library(tidyverse); library(RColorBrewer); library(deSolve)

# Function: Diffusion for Glucose and Microbes
RImodel <- function(time, state, parms, distance, nr0) {
  with(as.list(parms), {
    nr = nr0 #renaming this argument, throws error of 'nr' is used as argument name in ode.1D
    dr = distance/nr
    r <- seq(dr/2,distance - dr/2,dr) #all the midpoints/observations
    
    GLUCOSE <- state[1:N]
    MICROBE <- state[(N+1):(2*N)]
    
    #Glucose Diffusion Equation
    d1 = GLUCOSE[1:(nr-2)] - 2*GLUCOSE[2:(nr-1)] + GLUCOSE[3:nr]
    d2 = GLUCOSE[3:nr] - GLUCOSE[1:(nr-2)]
    dGLUCOSE = (d1/dr**2 + d2/2/dr/r[2:(nr-1)]) * D
    
    
    #Boundary Conditions
    b1 <- D/(dr^2) * (GLUCOSE[2] - GLUCOSE[1]) + D/(2*r[1]*dr) * (GLUCOSE[2] - GLUCOSE[1])
    b2 <- D/(dr^2) * (GLUCOSE[nr-1] - GLUCOSE[nr]) + D/(2*r[nr]*dr) * (GLUCOSE[nr-1] - GLUCOSE[nr])
    
    #Derivative vector 
    dGLUCOSE = c(b1, dGLUCOSE, b2)
    
    # Biology: Glucose and Microbes
    GlucoseUptake <- uptakeRate * GLUCOSE * MICROBE 
    MicrobeGrowth <- GlucoseUptake * CUE 
    MicrobeMortality <- mortalityRate * MICROBE
    depolymerization <-Vmaxprime * Microbe/(Microbe + km)
    
    # Rate of change = Flux gradient + Biology
    dGLUCOSE <-  dGLUCOSE - GlucoseUptake + MicrobeMortality + SOC_mineralization
    dMICROBE <- MicrobeGrowth - MicrobeMortality
    
    
    return(list(c(dGLUCOSE, dMICROBE)))
  })
}
########################################################################################################
#Time
dt = 3600/4 #time step (seconds)
days = 300#number of days to run the simulation
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
##########################################################################################################

#SOC Mineralization off 

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
           CUE = 0.5,
           SOC_mineralization = 0)


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

# account for the area 
# Calculate the total mass of glucose and microbes over time
area_matrix <- matrix(rep(area,length(time)),nrow=length(time),ncol=N,byrow=TRUE) #creates matrix with area repeated for each time

total_glucose_mass <-  rowSums(GLUCOSE_result * area_matrix) #get mass by multiplying concentration with area
total_microbe_mass <- rowSums(MICROBE_result * area_matrix)

#Area for each ring
area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi

#total area

total_area <-R^2 * pi 

# Aveg. Glucose & Microbes concentration over time
avg_glucose_concentration <- total_glucose_mass / total_area # calc. average concentrations over time 
avg_microbe_concentration <- total_microbe_mass / total_area

#time_vector <- time
plot( time/86400, avg_glucose_concentration, ylim = c(0,0.01),
      lines(time/86400, avg_microbe_concentration, col = "green"),xlab = "Days", ylab = "Concentration", 
      main = "Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "green"), lty = 1)

# Aveg. inner most ring Glucose & Microbes conc. across time 
inner_ring_glucose_concentration <- GLUCOSE_result[, 1] #Extract concentrations for the innermost ring
inner_ring_microbe_concentration <- MICROBE_result[, 1] 

plot(time/86400, inner_ring_glucose_concentration,
     lines(time/86400, inner_ring_microbe_concentration, col = "darkturquoise"),
     ylim = c(min(inner_ring_glucose_concentration, inner_ring_microbe_concentration), 
              max(inner_ring_glucose_concentration, inner_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Inner Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "darkturquoise"), lty = 1)

# Aveg. outer most ring Glucose & Microbes conc. across time 
outer_ring_glucose_concentration <- GLUCOSE_result[, N] 
outer_ring_microbe_concentration <- MICROBE_result[, N] 

plot(time/86400, outer_ring_glucose_concentration,
     lines(time/86400, outer_ring_microbe_concentration, col = "red"),
     ylim = c(min(outer_ring_glucose_concentration, outer_ring_microbe_concentration), 
              max(outer_ring_glucose_concentration, outer_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Outer Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "red"), lty = 1)

# Aveg. middle most ring Glucose & Microbes conc. across time 
middle_ring_glucose_concentration <- GLUCOSE_result [, (N + 1) %/% 2] 
middle_ring_microbe_concentration <- MICROBE_result [, (N + 1) %/% 2] 

plot(time/86400, middle_ring_glucose_concentration,
     lines(time/86400, middle_ring_microbe_concentration, col = "purple"),
     ylim = c(min(middle_ring_glucose_concentration, middle_ring_microbe_concentration), 
              max(middle_ring_glucose_concentration, middle_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Middle Ring Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "purple"), lty = 1)

# Distance of Microbe and Glucose at the end of the simulation at every location 
final_glucose_concentration <- GLUCOSE_result[nrow(GLUCOSE_result), ]
final_microbe_concentration <- MICROBE_result[nrow(MICROBE_result), ]

plot(r, final_glucose_concentration,
     xlab = "Location (r)", ylab = "Final Glucose Concentration", 
     main = "Final Glucose Concentration at the End of Simulation")

plot(r, final_microbe_concentration,
     xlab = "Location (r)", ylab = "Final Microbe Concentration", 
     main = "Final Microbe Concentration at the End of Simulation")

# Plot the evolution of microbial biomass at the initial location over time
plot(time/86400, MICROBE_result[, initial_microbe_location],
     xlab = "Days", ylab = "Microbial Biomass",
     main = "Microbial Biomass at Initial Location Over Time")


#Don't have any initial glucose & see the base release. ( Microbes bigger than 0)

# Model parameters
R <- 20
N <- 100
dr <- R/N
r <- seq(dr/2, by = dr, len = N)
ri <- seq(0, by = dr, len = N+1)
dri <- dr
parms <- c(D = 0.5 * 10^-4,
           uptakeRate = 0,
           mortalityRate = 0.2/86400,
           CUE = 0.5,
           SOC_mineralization = SOC_mineralization)

# Initial conditions
Glucose = rep(0, N)
#Glucose[1] <- C_add <- SOC * 0.01 
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

# account for the area 
# Calculate the total mass of glucose and microbes over time
area_matrix <- matrix(rep(area,length(time)),nrow=length(time),ncol=N,byrow=TRUE) #creates matrix with area repeated for each time

total_glucose_mass <-  rowSums(GLUCOSE_result * area_matrix) #get mass by multiplying concentration with area
total_microbe_mass <- rowSums(MICROBE_result * area_matrix)

#Area for each ring
area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi

#total area

total_area <-R^2 * pi 

# Aveg. Glucose & Microbes concentration over time
avg_glucose_concentration <- total_glucose_mass / total_area # calc. average concentrations over time 
avg_microbe_concentration <- total_microbe_mass / total_area

#time_vector <- time
plot( time/86400, avg_glucose_concentration, ylim = c(0,0.01),
      lines(time/86400, avg_microbe_concentration, col = "green"),xlab = "Days", ylab = "Concentration", 
      main = "Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "green"), lty = 1)

# Aveg. inner most ring Glucose & Microbes conc. across time 
inner_ring_glucose_concentration <- GLUCOSE_result[, 1] #Extract concentrations for the innermost ring
inner_ring_microbe_concentration <- MICROBE_result[, 1] 

plot(time/86400, inner_ring_glucose_concentration,
     lines(time/86400, inner_ring_microbe_concentration, col = "darkturquoise"),
     ylim = c(min(inner_ring_glucose_concentration, inner_ring_microbe_concentration), 
              max(inner_ring_glucose_concentration, inner_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Inner Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "darkturquoise"), lty = 1)

# Aveg. outer most ring Glucose & Microbes conc. across time 
outer_ring_glucose_concentration <- GLUCOSE_result[, N] 
outer_ring_microbe_concentration <- MICROBE_result[, N] 

plot(time/86400, outer_ring_glucose_concentration,
     lines(time/86400, outer_ring_microbe_concentration, col = "red"),
     ylim = c(min(outer_ring_glucose_concentration, outer_ring_microbe_concentration), 
              max(outer_ring_glucose_concentration, outer_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Outer Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "red"), lty = 1)

# Aveg. middle most ring Glucose & Microbes conc. across time 
middle_ring_glucose_concentration <- GLUCOSE_result [, (N + 1) %/% 2] 
middle_ring_microbe_concentration <- MICROBE_result [, (N + 1) %/% 2] 

plot(time/86400, middle_ring_glucose_concentration,
     lines(time/86400, middle_ring_microbe_concentration, col = "purple"),
     ylim = c(min(middle_ring_glucose_concentration, middle_ring_microbe_concentration), 
              max(middle_ring_glucose_concentration, middle_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Middle Ring Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "purple"), lty = 1)

# Distance of Microbe and Glucose at the end of the simulation at every location 
final_glucose_concentration <- GLUCOSE_result[nrow(GLUCOSE_result), ]
final_microbe_concentration <- MICROBE_result[nrow(MICROBE_result), ]

plot(r, final_glucose_concentration,
     xlab = "Location (r)", ylab = "Final Glucose Concentration", 
     main = "Final Glucose Concentration at the End of Simulation")

plot(r, final_microbe_concentration,
     xlab = "Location (r)", ylab = "Final Microbe Concentration", 
     main = "Final Microbe Concentration at the End of Simulation")

# Plot the evolution of microbial biomass at the initial location over time
plot(time/86400, MICROBE_result[, initial_microbe_location],
     xlab = "Days", ylab = "Microbial Biomass",
     main = "Microbial Biomass at Initial Location Over Time")


#Turn diffusion off & have an initial glucose concentration at the root surface with no base mineralization

# Model parameters
R <- 20
N <- 100
dr <- R/N
r <- seq(dr/2, by = dr, len = N)
ri <- seq(0, by = dr, len = N+1)
dri <- dr
parms <- c(D = 0,
           uptakeRate = 0.2/86400,
           mortalityRate = 0.2/86400,
           CUE = 0.5,
           SOC_mineralization = 0)

initial_glucose = initial_glucose

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

# account for the area 
# Calculate the total mass of glucose and microbes over time
area_matrix <- matrix(rep(area,length(time)),nrow=length(time),ncol=N,byrow=TRUE) #creates matrix with area repeated for each time

total_glucose_mass <-  rowSums(GLUCOSE_result * area_matrix) #get mass by multiplying concentration with area
total_microbe_mass <- rowSums(MICROBE_result * area_matrix)

#Area for each ring
area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi

#total area

total_area <-R^2 * pi 

# Aveg. Glucose & Microbes concentration over time
avg_glucose_concentration <- total_glucose_mass / total_area # calc. average concentrations over time 
avg_microbe_concentration <- total_microbe_mass / total_area

#time_vector <- time
plot( time/86400, avg_glucose_concentration, ylim = c(0,0.01),
      lines(time/86400, avg_microbe_concentration, col = "green"),xlab = "Days", ylab = "Concentration", 
      main = "Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "green"), lty = 1)

# Aveg. inner most ring Glucose & Microbes conc. across time 
inner_ring_glucose_concentration <- GLUCOSE_result[, 1] #Extract concentrations for the innermost ring
inner_ring_microbe_concentration <- MICROBE_result[, 1] 

plot(time/86400, inner_ring_glucose_concentration,
     lines(time/86400, inner_ring_microbe_concentration, col = "darkturquoise"),
     ylim = c(min(inner_ring_glucose_concentration, inner_ring_microbe_concentration), 
              max(inner_ring_glucose_concentration, inner_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Inner Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "darkturquoise"), lty = 1)

# Aveg. outer most ring Glucose & Microbes conc. across time 
outer_ring_glucose_concentration <- GLUCOSE_result[, N] 
outer_ring_microbe_concentration <- MICROBE_result[, N] 

plot(time/86400, outer_ring_glucose_concentration,
     lines(time/86400, outer_ring_microbe_concentration, col = "red"),
     ylim = c(min(outer_ring_glucose_concentration, outer_ring_microbe_concentration), 
              max(outer_ring_glucose_concentration, outer_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Outer Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "red"), lty = 1)

# Aveg. middle most ring Glucose & Microbes conc. across time 
middle_ring_glucose_concentration <- GLUCOSE_result [, (N + 1) %/% 2] 
middle_ring_microbe_concentration <- MICROBE_result [, (N + 1) %/% 2] 

plot(time/86400, middle_ring_glucose_concentration,
     lines(time/86400, middle_ring_microbe_concentration, col = "purple"),
     ylim = c(min(middle_ring_glucose_concentration, middle_ring_microbe_concentration), 
              max(middle_ring_glucose_concentration, middle_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Middle Ring Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "purple"), lty = 1)

# Distance of Microbe and Glucose at the end of the simulation at every location 
final_glucose_concentration <- GLUCOSE_result[nrow(GLUCOSE_result), ]
final_microbe_concentration <- MICROBE_result[nrow(MICROBE_result), ]

plot(r, final_glucose_concentration,
     xlab = "Location (r)", ylab = "Final Glucose Concentration", 
     main = "Final Glucose Concentration at the End of Simulation")

plot(r, final_microbe_concentration,
     xlab = "Location (r)", ylab = "Final Microbe Concentration", 
     main = "Final Microbe Concentration at the End of Simulation")

# Plot the evolution of microbial biomass at the initial location over time
plot(time/86400, MICROBE_result[, initial_microbe_location],
     xlab = "Days", ylab = "Microbial Biomass",
     main = "Microbial Biomass at Initial Location Over Time")

# Run with just diffusion to see mass balance
# Model parameters
R <- 20
N <- 100
dr <- R/N
r <- seq(dr/2, by = dr, len = N)
ri <- seq(0, by = dr, len = N+1)
dri <- dr
parms <- c(D = 0.5 * 10^-4,
           uptakeRate = 0,
           mortalityRate = 0,
           CUE = 0,
           SOC_mineralization = 0)

#Area for each ring
area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi

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


# Calculate the total mass of glucose and microbes over time
area_matrix <- matrix(rep(area,length(time)),nrow=length(time),ncol=N,byrow=TRUE) #creates matrix with area repeated for each time

total_glucose_mass <-  rowSums(GLUCOSE_result * area_matrix) #get mass by multiplying concentration with area
total_microbe_mass <- rowSums(MICROBE_result * area_matrix)


plot(time, total_glucose_mass)

#Area for each ring
area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi

#total area

total_area <-R^2 * pi 

# Aveg. Glucose & Microbes concentration over time
avg_glucose_concentration <- total_glucose_mass / total_area # calc. average concentrations over time 
avg_microbe_concentration <- total_microbe_mass / total_area

# Aveg. Glucose & Microbes concentration over time
avg_glucose_concentration <- total_glucose_mass / total_area # calc. average concentrations over time 
avg_microbe_concentration <- total_microbe_mass / total_area

#time_vector <- time
plot( time/86400, avg_glucose_concentration, ylim = c(0,0.01),
      lines(time/86400, avg_microbe_concentration, col = "green"),xlab = "Days", ylab = "Concentration", 
      main = "Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "green"), lty = 1)

# Aveg. inner most ring Glucose & Microbes conc. across time 
inner_ring_glucose_concentration <- GLUCOSE_result[, 1] #Extract concentrations for the innermost ring
inner_ring_microbe_concentration <- MICROBE_result[, 1] 

plot(time/86400, inner_ring_glucose_concentration,
     lines(time/86400, inner_ring_microbe_concentration, col = "darkturquoise"),
     ylim = c(min(inner_ring_glucose_concentration, inner_ring_microbe_concentration), 
              max(inner_ring_glucose_concentration, inner_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Inner Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "darkturquoise"), lty = 1)

# Aveg. outer most ring Glucose & Microbes conc. across time 
outer_ring_glucose_concentration <- GLUCOSE_result[, N] 
outer_ring_microbe_concentration <- MICROBE_result[, N] 

plot(time/86400, outer_ring_glucose_concentration,
     lines(time/86400, outer_ring_microbe_concentration, col = "red"),
     ylim = c(min(outer_ring_glucose_concentration, outer_ring_microbe_concentration), 
              max(outer_ring_glucose_concentration, outer_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Outer Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "red"), lty = 1)

# Aveg. middle most ring Glucose & Microbes conc. across time 
middle_ring_glucose_concentration <- GLUCOSE_result [, (N + 1) %/% 2] 
middle_ring_microbe_concentration <- MICROBE_result [, (N + 1) %/% 2] 

plot(time/86400, middle_ring_glucose_concentration,
     lines(time/86400, middle_ring_microbe_concentration, col = "purple"),
     ylim = c(min(middle_ring_glucose_concentration, middle_ring_microbe_concentration), 
              max(middle_ring_glucose_concentration, middle_ring_microbe_concentration)),
     xlab = "Days", ylab = "Concentration", 
     main = "Middle Ring Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "purple"), lty = 1)

# Distance of Microbe and Glucose at the end of the simulation at every location 
final_glucose_concentration <- GLUCOSE_result[nrow(GLUCOSE_result), ]
final_microbe_concentration <- MICROBE_result[nrow(MICROBE_result), ]

plot(r, final_glucose_concentration,
     xlab = "Location (r)", ylab = "Final Glucose Concentration", 
     main = "Final Glucose Concentration at the End of Simulation")

plot(r, final_microbe_concentration,
     xlab = "Location (r)", ylab = "Final Microbe Concentration", 
     main = "Final Microbe Concentration at the End of Simulation")

# Plot the evolution of microbial biomass at the initial location over time
plot(time/86400, MICROBE_result[, initial_microbe_location],
     xlab = "Days", ylab = "Microbial Biomass",
     main = "Microbial Biomass at Initial Location Over Time")

#----------------------------------------------------------------------------------

SOC_concentration <- bulk_density * 0.02 #mg
SOC_mineralization_rate <- 0.08/365/86400  #second -1 (replace with your actual value = .00002)
SOC_mineralization <- SOC_concentration * SOC_mineralization_rate #mg * mm-3 * s-1

# Model parameters (no uptake, no mortality, no CUE, no SOC mineralization)
R <- 20
N <- 100
dr <- R/N
r <- seq(dr/2, by = dr, len = N)
ri <- seq(0, by = dr, len = N+1)
parms <- c(D = 0.5 * 10^-4,
           uptakeRate = 300/86400,
           mortalityRate = 0.033/86400,
           CUE = 0.5,
           SOC_mineralization = SOC_mineralization)

# Initial conditions
C_add <- SOC * 0.01 #Add 1% of SOC in terms of glucose carbon (mg)
Glucose = rep(C_add/1000, N) #same everywhere

# Set the initial location of the microbial biomass (constant)
Microbes <- rep(C_add/100, N)

state <- c(Glucose, Microbes)

# Solve the model using the default banded reformulation method
model <- RImodel(0, state, parms, distance, N)
result <- ode.1D(y = state, times = time, func = RImodel, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time <- result[, 1]

# Extracting results
GLUCOSE_result <- result[, 2:(N + 1)]
MICROBE_result <- result[, (N + 2):(2*N + 1)]

# Calculate the total mass of glucose and microbes over time 
area_matrix <- matrix(rep(area,length(time)),nrow=length(time),ncol=N,byrow=TRUE) #creates matrix with area repeated for each time

total_glucose_mass <-  rowSums(GLUCOSE_result * area_matrix) #get mass by multiplying concentration with area
total_microbe_mass <- rowSums(MICROBE_result * area_matrix)


#Area for each ring
area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi

#total area

total_area <-R^2 * pi 

# Aveg. Glucose & Microbes concentration over time
avg_glucose_concentration <- total_glucose_mass / total_area # calc. average concentrations over time 
avg_microbe_concentration <- total_microbe_mass / total_area

#time_vector <- time
plot( time/86400, avg_glucose_concentration, ylim = c(0,0.01),
lines(time/86400, avg_microbe_concentration, col = "green"),xlab = "Days", ylab = "Concentration", 
main = "Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "green"), lty = 1)

# Aveg. inner most ring Glucose & Microbes conc. across time 
inner_ring_glucose_concentration <- GLUCOSE_result[, 1] #Extract concentrations for the innermost ring
inner_ring_microbe_concentration <- MICROBE_result[, 1] 

plot(time/86400, inner_ring_glucose_concentration,
lines(time/86400, inner_ring_microbe_concentration, col = "darkturquoise"),
ylim = c(min(inner_ring_glucose_concentration, inner_ring_microbe_concentration), 
         max(inner_ring_glucose_concentration, inner_ring_microbe_concentration)),
xlab = "Days", ylab = "Concentration", 
main = "Inner Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "darkturquoise"), lty = 1)

# Aveg. outer most ring Glucose & Microbes conc. across time 
outer_ring_glucose_concentration <- GLUCOSE_result[, N] 
outer_ring_microbe_concentration <- MICROBE_result[, N] 

plot(time/86400, outer_ring_glucose_concentration,
lines(time/86400, outer_ring_microbe_concentration, col = "red"),
ylim = c(min(outer_ring_glucose_concentration, outer_ring_microbe_concentration), 
         max(outer_ring_glucose_concentration, outer_ring_microbe_concentration)),
xlab = "Days", ylab = "Concentration", 
main = "Outer Ring Glucose and Microbe Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "red"), lty = 1)

# Aveg. middle most ring Glucose & Microbes conc. across time 
middle_ring_glucose_concentration <- GLUCOSE_result [, (N + 1) %/% 2] 
middle_ring_microbe_concentration <- MICROBE_result [, (N + 1) %/% 2] 

plot(time/86400, middle_ring_glucose_concentration,
lines(time/86400, middle_ring_microbe_concentration, col = "purple"),
ylim = c(min(middle_ring_glucose_concentration, middle_ring_microbe_concentration), 
         max(middle_ring_glucose_concentration, middle_ring_microbe_concentration)),
xlab = "Days", ylab = "Concentration", 
main = "Middle Ring Concentration Over Time")

legend("topright", legend = c("Glucose Conc.", "Microbe Conc."), col = c("black", "purple"), lty = 1)

# Distance of Microbe and Glucose at the end of the simulation at every location 
final_glucose_concentration <- GLUCOSE_result[nrow(GLUCOSE_result), ]
final_microbe_concentration <- MICROBE_result[nrow(MICROBE_result), ]

plot(r, final_glucose_concentration,
xlab = "Location (r)", ylab = "Final Glucose Concentration", 
main = "Final Glucose Concentration at the End of Simulation")

plot(r, final_microbe_concentration,
xlab = "Location (r)", ylab = "Final Microbe Concentration", 
main = "Final Microbe Concentration at the End of Simulation")

# Plot the evolution of microbial biomass at the initial location over time
plot(time/86400, MICROBE_result[, initial_microbe_location],
xlab = "Days", ylab = "Microbial Biomass",
main = "Microbial Biomass at Initial Location Over Time")































