#Load in libraries ----
library(tidyverse); library(RColorBrewer); library(deSolve)

# Function: Diffusion for Glucose and Microbes
RImodel <- function(time, state, parms, distance, nr0) {
  with(as.list(parms), {
    nr = nr0 # renaming this argument, throws error of 'nr' is used as argument name in ode.1D
    dr = distance/nr
    r <- seq(dr/2,distance - dr/2,dr) # all the midpoints/observations
    
    #Extracting Glucose and Microbe concentrations from the 'state' vector
    GLUCOSE <- state[1:N]   
    MICROBE <- state[(N+1):(2*N)]
    
    #Glucose Diffusion Equation
    d1 = GLUCOSE[1:(nr-2)] - 2*GLUCOSE[2:(nr-1)] + GLUCOSE[3:nr] # Computing the first and second derivatives for the diffusion equation using finite differences
    d2 = GLUCOSE[3:nr] - GLUCOSE[1:(nr-2)]
    
    #Computing the derivative for glucose concentration considering diffusion (D)
    dGLUCOSE = (d1/dr**2 + d2/2/dr/r[2:(nr-1)]) * D
    
    #Boundary Conditions
    # Applying boundary conditions for the diffusion equation at the start and end points
    b1 <- D/(dr^2) * (GLUCOSE[2] - GLUCOSE[1]) + D/(2*r[1]*dr) * (GLUCOSE[2] - GLUCOSE[1]) 
    b2 <- D/(dr^2) * (GLUCOSE[nr-1] - GLUCOSE[nr]) + D/(2*r[nr]*dr) * (GLUCOSE[nr-1] - GLUCOSE[nr])
    
    #Derivative vector for Glucose Diffusion Equation
    dGLUCOSE = c(b1, dGLUCOSE, b2)  # Concatenating boundary conditions with the computed derivatives
    
    #Biology: Glucose and Microbes (Exchange)
    GlucoseUptake <- uptakeRate * GLUCOSE * MICROBE     # GlucoseUptake represents the rate at which microbes take up glucose from their environment.It is calculated by multiplying the uptakeRate with the current concentrations of GLUCOSE and MICROBE.
    MicrobeGrowth <- GlucoseUptake * CUE                # MicrobeGrowth represents the growth of microbes based on glucose uptake and Carbon Use Efficiency (CUE).It quantifies the increase in microbe population due to the availability of glucose, considering the efficiency of resource utilization (CUE).
    MicrobeMortality <- mortalityRate * MICROBE         # MicrobeMortality represents the rate at which microbes die off. It is calculated by multiplying the mortalityRate with the current concentration of MICROBE.
    depolymerization <-Vmaxprime * MICROBE/(MICROBE+km) # Depolymerization represents the breakdown of complex molecules into simpler units by microbes.Vmaxprime is the maximum rate of depolymerization.
                                                        # MICROBE is the microbe concentration, and km is the Michaelis constant. The rate of depolymerization increases with higher MICROBE concentration but levels off at higher concentrations due to the Michaelis constant.
    
    #Rate of change = Flux gradient + Biology 
    dGLUCOSE <-  dGLUCOSE - GlucoseUptake + MicrobeMortality + depolymerization # The change in glucose concentration over time is influenced by: Diffusion processes causing glucose to move within the system,Glucose uptake by microbes-reducing glucose concentration
                                                                                # Microbe mortality-leading to the release of glucose back into the environment & Microbial depolymerization, breaking down complex molecules into glucose 
    dMICROBE <- MicrobeGrowth - MicrobeMortality                                # The change in microbe population over time is influenced by: Microbial growth due to glucose uptake & Microbe mortality, representing the death of microbes 
    
    
    return(list(c(dGLUCOSE, dMICROBE))) # Returning the derivative vectors for glucose diffusion and microbe dynamics.
  })
}

# Time
dt = 3600/4 #time step (seconds)
days = 14  #number of days to run the simulation
duration = days*86400 #days to seconds
nt = duration/dt +1 #number of observations
time = seq(0,duration,dt) #seconds
 
# Space (1-D) 
r0 = 0.03 #where root meets soil (mm)
distance = 30 #distance away from center (mm) 

# Parameters
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
# Model application

# Model parameters
R <- 20     #Radius of the system
N <- 100    #Number of spatial grid points
dr <- R/N   #Spatial grid size
r <- seq(dr/2, by = dr, len = N)
ri <- seq(0, by = dr, len = N+1) #Inner radius of each grid point
dri <- dr   #Grid spacing
parms <- c(D = 0.5e-4,             #Diffusion coefficient
           uptakeRate = 2.4e-11,  #Initial uptake rate of glucose by microbes
           mortalityRate = 1.05e-9,   #Microbial mortality rate
           CUE = 0.5,                   #Carbon use efficiency
           Vmaxprime = 3.7e-7,         #Maximum microbial growth rate
           km = 88)                         #Michaelis-Menten constant for glucose uptake


# Initial conditions

#The initial glucose concentration at the first grid point is calculated as 1% of SOC (Soil Organic Carbon) concentration
Glucose = rep(0, N) #Initialize an array of length N (number of grid points) with initial glucose concentration set to 0
Glucose[1] <- C_add <- SOC * 0.01 #Set the initial glucose concentration at the first grid point
Microbes = rep(C_add/100, N) #The initial microbial concentration is set to a fraction (1/100) of the initial glucose concentration (C_add).
                             #Initialize an array of length N with initial microbial concentration set to a fraction of C_add

#'state' is a concatenated vector containing the initial concentrations of glucose followed by microbial concentrations.
state = c(Glucose,Microbes) #Combine the arrays of initial glucose and microbial concentrations into a single state vector

# Time points for output

# Solve the model using the default banded reformulation method
model <- RImodel(0, state, parms, distance, N)
result <- ode.1D(y = state, times = time, func = RImodel, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time = result[,1]

# Extracting results
GLUCOSE_result <- result[, 2:(N + 1)]
MICROBE_result <- result[, (N + 2):(2*N + 1)]

# Account for the area 
# Calculate the total mass of glucose and microbes over time
area_matrix <- matrix(rep(area,length(time)),nrow=length(time),ncol=N,byrow=TRUE) #creates matrix with area repeated for each time

total_glucose_mass <-  rowSums(GLUCOSE_result * area_matrix) #get mass by multiplying concentration with area
total_microbe_mass <- rowSums(MICROBE_result * area_matrix)

#Area for each ring
area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi

#total area
total_area <-R^2 * pi 

# Plotting the total mass of glucose over time
plot(time/86400, total_glucose_mass, 
     type = "l", 
     xlab = "Time", 
     ylab = "Total Glucose Mass (g)",
     main = "Total Glucose Mass Over Time",

col = "orange",                              # Set line color 
ylim = c(0, max(total_glucose_mass) * 1.2))  # Set y-axis limit to 20% more than the maximum value

# Adding labels and annotations
abline(h = max(total_glucose_mass), col = "red", lty = 2)  # Add a horizontal dashed line at the maximum glucose mass
text(max(time/86400), max(total_glucose_mass), 
     paste("Max Glucose Mass:", round(max(total_glucose_mass), 2), "g"), 
     pos = 4, col = "red")  # Annotate the maximum glucose mass value

# Plotting the total mass of microbes over time
plot(time/86400, total_microbe_mass, 
     type = "l", 
     xlab = "Time", 
     ylab = "Total Microbe Mass", 
     main = "Total Microbe Mass Over Time",

col = "green",     # Set line color to green
ylim = c(0, max(total_microbe_mass) * 1.2))  # Set y-axis limit to 20% more than the maximum value

# Adding labels and annotations
abline(h = max(total_microbe_mass), col = "red", lty = 2)  # Add a horizontal dashed line at the maximum microbe mass
text(max(time/86400), max(total_microbe_mass), 
     paste("Max Microbe Mass:", round(max(total_microbe_mass), 2), "g"), 
     pos = 4, col = "red")  # Annotate the maximum microbe mass value

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

