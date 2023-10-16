#Model Expansion 

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
    
    # Rate of change = Flux gradient + Biology
    dGLUCOSE <-  dGLUCOSE - GlucoseUptake + MicrobeMortality + SOC_mineralization
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
           CUE = 0.5,
           SOC_mineralization = SOC_mineralization)


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

###########################################################################

#Add SOC Mineralization 
#Add Base Mineralization Supply 

# Parameters
bulk_density <- 1.01  # g/cm^3 = mg/mm^3
soil_mass <- bulk_density * area_volume  # mg
SOC <- soil_mass * 0.02  # Assumes SOC is 2% (mg)
SOC_concentration <- bulk_density * 0.02 #mg

#Linear relationship between SOC and mineralization rate
mineralization_rate_per_percent <- 0.001  # Example factor (adjust as needed)
SOC_mineralization_rate <- SOC * mineralization_rate_per_percent

# Print the calculated SOC mineralization rate
cat("SOC Mineralization Rate:", SOC_mineralization_rate, "mg/day\n")

# Internal mineralization/SOC
SOC_mineralization_rate <- 0.05/365/86400  #second -1 (replace with your actual value = .00002)
SOC_mineralization <- SOC_concentration * SOC_mineralization_rate #mg * mm-3 * s-1

 # Calculate total mass of glucose and microbes for each timestep
 glucose_total_mass <- colSums(GLUCOSE_result) * dri
 microbes_total_mass <- colSums(MICROBE_result) * dri
 
 # Convert the total masses to milligrams (assuming your model uses mg/mm as units)
 glucose_total_mass <- glucose_total_mass * area_volume
 microbes_total_mass <- microbes_total_mass * area_volume
 
 #----------------------#
      #Mass Balance
 #----------------------#
 
 #Area for each ring
 area = (ri[2:(N+1)]^2 - ri[1:(N)]^2) * pi
 
 # Calculate the initial mass of glucose & microbes
 initial_glucose_mass <- sum(Glucose * area) 
 initial_microbes_mass <- sum(Microbes * area) 
 
 # Variables to track cumulative mass changes
cum_glucose_mass_change <- numeric(length(time))
cum_microbes_mass_change <- numeric(length(time))

 # Calculate cumulative mass changes for each timestep 
for (i in 2:length(time)) {

  glucose_mass_change <- sum((GLUCOSE_result[i,] * area - GLUCOSE_result[i-1,] * area) )
  cum_glucose_mass_change[i] <- glucose_mass_change
  
  microbes_mass_change <- sum((MICROBE_result[i,] * area - MICROBE_result[i-1,] * area))
  cum_microbes_mass_change[i] <- microbes_mass_change
}


# Calculate mass balance 
final_glucose_mass <-initial_glucose_mass + sum(cum_glucose_mass_change)
final_microbes_mass <-initial_microbes_mass + sum(cum_microbes_mass_change)

# Print mass balance results
cat("Initial Glucose Mass:", initial_glucose_mass, "mg\n")
cat("Final Glucose Mass:", final_glucose_mass, "mg\n")
cat("Initial Microbes Mass:", initial_microbes_mass, "mg\n")
cat("Final Microbes Mass:", final_microbes_mass, "mg\n")
 
 #-------------------------------#
  #Total mass for each time step
 #-------------------------------#

 # Variables to store total mass for each time step
total_mass_glucose <- numeric(length(time))
total_mass_microbes <- numeric(length(time))

 # Calculate total mass for each time step
for (i in 1:length(time)) {
  total_mass_glucose[i] <- sum(GLUCOSE_result[i,] * area) 
  total_mass_microbes[i] <- sum(MICROBE_result[i,] * area) 
}

 # Print total mass for each time step
for (i in 1:length(time)) {
  cat("Time:", time[i], "s\n")
  cat("Total Glucose Mass:", total_mass_glucose[i], "mg\n")
  cat("Total Microbes Mass:", total_mass_microbes[i], "mg\n")
  cat("--------------------\n")
  
}

 #-------------------------------------------------------------------#
  #SOC_Mineralization comparison between 0 & the value given already 
 #-------------------------------------------------------------------#
 
# Calculate total SOC Mineralization for each time step
SOC_mineralization_total <- colSums(GLUCOSE_result * SOC_mineralization_rate) * area

# Calculate SOC Mineralization for each time step when SOC is 0 
SOC_mineralization_zero <- rep(0, length(time)) # Assuming no mineralization when SOC = 0

# Calculate the cumulative SOC mineralization for the entire simulation
cum_mineralization_total <- sum(SOC_mineralization_total)
cum_mineralization_zero <- sum(SOC_mineralization_zero)

# Compare cumulative SOC Mineralization 
if (cum_mineralization_total > cum_mineralization_zero) {
  cat("Cumulative SOC Mineralization is greater than 0\n")
} else if (cum_mineralization_total < cum_mineralization_zero) {
  cat("Cumulative SOC Mineralization is less than 0\n")
} else {
  cat("Cumulative SOC Mineralization is equal to 0\n")
}
#calculate the cumulative SOC mineralization over the entire simulation period for both cases (given value and when SOC is 0). Then, you compare the cumulative values to determine whether mineralization was greater, less, or equal to 0.
#_____________________________

# Calculate total SOC Mineralization for each time step
SOC_mineralization_total <- GLUCOSE_result * SOC_mineralization_rate * area

# Calculate SOC Mineralization for each time step when SOC is 0
SOC_mineralization_zero <- matrix(0, nrow = length(time), ncol = N)  # Assuming no mineralization when SOC = 0

# Variables to store comparative results and numeric values for each timestep
comparative_results <- character(length(time))
numeric_values <- numeric(length(time))

# Loop through each timestep
for (i in 1:length(time)) {
  numeric_values[i] <- sum(SOC_mineralization_total[i,])
  
  if (numeric_values[i] > sum(SOC_mineralization_zero[i,])) {
    comparative_results[i] <- "greater than 0"
  } else if (numeric_values[i] < sum(SOC_mineralization_zero[i,])) {
    comparative_results[i] <- "less than 0"
  } else {
    comparative_results[i] <- "equal to 0"
  }
}

# Print comparative results and numeric values for each timestep
for (i in 1:length(time)) {
  cat("Time:", time[i], "s\n")
  cat("Total SOC Mineralization:", numeric_values[i], "\n")
  cat("Comparative Result:", comparative_results[i], "\n")
  cat("--------------------\n")
}

###########################################################################

 #Calculate change in microbial populations 
microbes_change <- numeric(length(time))
microbes_change[1] <- 0

for (i in 2:length(time)) {
  microbes_change[i] <- sum(abs(MICROBE_result[i,] - MICROBE_result[i - 1,]))
}

# Determine Equilibrium for microbes
equilibrium_time <- NULL
equilibrium_threshold <- 1e-6  # Set your desired threshold for equilibrium

for (i in 2:length(time)) {
  if (microbes_change[i] < equilibrium_threshold) {
    equilibrium_time <- time[i]
    break  # Stop loop when equilibrium is reached
  }
}

if (!is.null(equilibrium_time)) {
  cat("Microbes reached equilibrium at time:", equilibrium_time, "s\n")
} else {
  cat("Microbes did not reach equilibrium within the simulation period\n")
}

 # Run model until microbes reach equilibrium !!!!!!!!!!!!!!!!!!!!!!!!!!!
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
    
    # Rate of change = Flux gradient + Biology
    dGLUCOSE <-  dGLUCOSE - GlucoseUptake + MicrobeMortality + SOC_mineralization
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
           CUE = 0.5,
           SOC_mineralization = SOC_mineralization)


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
equilibrium_reached = FALSE

# Loop for reaching equilibrium
while (!equilibrium_reached) {
  # Update the model using the current state
  model <- RImodel(time[length(time)], state, parms, distance, N)
  
  # Calculate the microbial population change
  microbes_change <- sum(abs(MICROBE_result[length(time),] - MICROBE_result[length(time) - 1,]))
  
  # Check for equilibrium
  if (microbes_change < equilibrium_threshold) {
    equilibrium_reached <- TRUE
    equilibrium_time <- time[length(time)]
  } else {
    # Extend the simulation time
    extended_time <- seq(time[length(time)], time[length(time)] + 10 * 86400, by=86400)
    Glucose = GLUCOSE_result[length(time),]
    Microbes =MICROBE_result[length(time),]
    state = c(Glucose,Microbes)
    extended_result <- ode.1D(y = state, times = extended_time, func = RImodel, parms = parms, distance = distance, nr0 = N,
                              nspec = 2, dimens = N)
    
    # Append extended results to existing results
    time <- c(time, extended_time)
    GLUCOSE_result <- rbind(GLUCOSE_result, extended_result[, 1:N])
    MICROBE_result <- rbind(MICROBE_result, extended_result[, (N + 1):(2*N)])
  }
}

# Print equilibrium results
cat("Microbes reached equilibrium at time:", equilibrium_time, "s\n")

 # Create 2nd Run & have this as initial conditions 

# Identify the time point where microbes reached equilibrium
equilibrium_index <- which(time == equilibrium_time)

# Extract the state at equilibrium
equilibrium_state <- c(GLUCOSE_result[equilibrium_index,], MICROBE_result[equilibrium_index,])

# Set the equilibrium state as initial conditions for the second simulation
state_modified <- equilibrium_state

# Run the second simulation with modified initial conditions
result_modified <- ode.1D(y = state_modified, times = time, func = RImodel, parms = parms, distance = R, nr0 = N, dimens = N,
                          nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time_modified <- result_modified[, 1]
GLUCOSE_result_modified <- result_modified[, 2:(N + 1)]
MICROBE_result_modified <- result_modified[, (N + 2):(2*N + 1)]
##########################################################################################

#Run with just diffusion to see mass balance!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # Load required libraries
library(deSolve)

# Function: Diffusion for Glucose
DiffusionModel <- function(time, state, parms, distance, nr0) {
  with(as.list(parms), {
    nr = nr0
    dr = distance/nr
    r <- seq(dr/2, distance - dr/2, dr)
    
    GLUCOSE <- state[1:nr]
    
    # Glucose Diffusion Equation
    d1 = GLUCOSE[1:(nr-2)] - 2*GLUCOSE[2:(nr-1)] + GLUCOSE[3:nr]
    d2 = GLUCOSE[3:nr] - GLUCOSE[1:(nr-2)]
    dGLUCOSE = (d1/(dr^2) + d2/(2*dr*r[2:(nr-1)])) * D
    
    # Boundary Conditions
    b1 <- D/(dr^2) * (GLUCOSE[2] - GLUCOSE[1]) + D/(2*r[1]*dr) * (GLUCOSE[2] - GLUCOSE[1])
    b2 <- D/(dr^2) * (GLUCOSE[nr-1] - GLUCOSE[nr]) + D/(2*r[nr]*dr) * (GLUCOSE[nr-1] - GLUCOSE[nr])
    
    # Derivative vector 
    dGLUCOSE <- c(b1, dGLUCOSE, b2)
    
    return(list(dGLUCOSE))
  })
}

# Time
dt <- 3600/4
days <- 10
duration <- days * 86400
nt <- duration/dt + 1
time <- seq(0, duration, dt)

# Space (1-D)
r0 <- 0.03
distance <- 30
R <- distance
N <- 100
dr <- R/N
parms <- list(D = 0.5 * 10^-4)

# Initial condition for glucose concentration
Glucose <- rep(0, N)
Glucose[1] <- 0.01 * 0.02 * area_volume  # Initial glucose mass (assuming SOC is 2% and 1% of SOC is added)

state <- Glucose

# Solve the diffusion model
result <- ode(y = state, times = time, func = DiffusionModel, parms = parms, distance = R, nr0 = N)

# Calculate total mass of glucose for each timestep

area_volume = (ri[2:(N+1)]^2 - ri[1:N]^2)*pi
total_mass_glucose <- numeric(nt)
nrow_results = length(result[,1])

for (i in 1:nrow_results){
  total_mass_glucose[i] = sum(result[i,-1]*area_volume)
}

mass_change <- diff(total_mass_glucose) #change for e/a time step

# Calculate cumulative mass change for each timestep
cumulative_mass_change <- cumsum(mass_change)

plot(time[1:(nrow_results - 1)],mass_change)


#Turn diffusion off & have an initial glucose concentration at the root surface with no base mineralization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Load required libraries
library(deSolve)

# Function: RImodel with no diffusion and no base mineralization
RImodel <- function(time, state, parms, distance, nr0) {
  with(as.list(parms), {
    nr = nr0
    GLUCOSE <- state[1:N]
    MICROBE <- state[(N+1):(2*N)]
    
    dGLUCOSE <- rep(0, N)
    dMICROBE <- rep(0, N)
    
    GlucoseUptake <- uptakeRate * GLUCOSE * MICROBE 
    MicrobeGrowth <- GlucoseUptake * CUE 
    MicrobeMortality <- mortalityRate * MICROBE
    
    dGLUCOSE <- dGLUCOSE - GlucoseUptake + MicrobeMortality
    dMICROBE <- MicrobeGrowth - MicrobeMortality
    
    return(list(c(dGLUCOSE, dMICROBE)))
  })
}

#------------------------------------------------------------------------------------------------------------------------
# Time
dt <- 3600/4  # time step (seconds)
days <- 10    # number of days to run the simulation
duration <- days * 86400  # days to seconds
nt <- duration/dt + 1     # number of observations
time <- seq(0, duration, dt)  # seconds

# Space (1-D)
r0 <- 0.03    # where root meets soil (mm)
distance <- 30  # distance away from center (mm)

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
#--------------------------------------------------------------------------------------------------------------------------
# Model parameters and initial conditions
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

# Initial conditions
initial_glucose_concentration <- 100  # Adjust as needed
initial_microbe_concentration <- 0.1   # Adjust as needed

Glucose <- rep(0, N)
Glucose[1] <- initial_glucose_concentration
Microbes <- rep(initial_microbe_concentration, N)

state <- c(Glucose, Microbes)

# Solve the model
model <- RImodel(0, state, parms, distance, N)
result <- ode.1D(y = state, times = time, func = RImodel, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time = result[, 1]
GLUCOSE_result <- result[, 2:(N + 1)]
MICROBE_result <- result[, (N + 2):(2*N + 1)]


# Print the simulation results for glucose and microbe concentrations
print(result)


#Dont have any initial glucose & see the base release. ( Microbes bigger than 0)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

# Load required libraries
library(deSolve)

# Function: RImodel with no initial glucose and base release
RImodel <- function(time, state, parms, distance, nr0) {
  with(as.list(parms), {
    nr = nr0
    GLUCOSE <- state[1:N]
    MICROBE <- state[(N+1):(2*N)]
    
    dGLUCOSE <- rep(0, N)
    dMICROBE <- rep(0, N)
    
    GlucoseUptake <- uptakeRate * GLUCOSE * MICROBE 
    MicrobeGrowth <- GlucoseUptake * CUE 
    MicrobeMortality <- mortalityRate * MICROBE
    
    dGLUCOSE <- dGLUCOSE - GlucoseUptake + MicrobeMortality + base_release
    dMICROBE <- MicrobeGrowth - MicrobeMortality
    
    return(list(c(dGLUCOSE, dMICROBE)))
  })
}

#------------------------------------------------------------------------------------------------------------------------
# Time
dt <- 3600/4  # time step (seconds)
days <- 10    # number of days to run the simulation
duration <- days * 86400  # days to seconds
nt <- duration/dt + 1     # number of observations
time <- seq(0, duration, dt)  # seconds

# Space (1-D)
r0 <- 0.03    # where root meets soil (mm)
distance <- 30  # distance away from center (mm)

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
#--------------------------------------------------------------------------------------------------------------------------

# Model parameters and initial conditions
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
           SOC_mineralization = SOC_mineralization)

# Other parameter values and calculations
base_release <- 0.1  # Adjust as needed

# Initial conditions
Glucose <- rep(0, N)
Microbes <- rep(0.1, N)  # Adjust initial microbe concentration if needed

state <- c(Glucose, Microbes)

# Solve the model
model <- RImodel(0, state, parms, distance, N)
result <- ode.1D(y = state, times = time, func = RImodel, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time <- result[, 1]
GLUCOSE_result <- result[, 2:(N + 1)]
MICROBE_result <- result[, (N + 2):(2*N + 1)]

# Print the simulation results for glucose and microbe concentrations
print(result)


#Add back diffusion & No Glucose!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(deSolve)

# Function: RImodel with diffusion and no initial glucose
RImodel <- function(time, state, parms, distance, nr0) {
  with(as.list(parms), {
    nr = nr0
    dr = distance/nr
    r <- seq(dr/2, distance - dr/2, dr)
    
    GLUCOSE <- state[1:N]
    MICROBE <- state[(N+1):(2*N)]
    
    dGLUCOSE <- rep(0, N)
    dMICROBE <- rep(0, N)
    
    d1 = GLUCOSE[1:(nr-2)] - 2*GLUCOSE[2:(nr-1)] + GLUCOSE[3:nr]
    d2 = GLUCOSE[3:nr] - GLUCOSE[1:(nr-2)]
    dGLUCOSE = (d1/dr**2 + d2/2/dr/r[2:(nr-1)]) * D
    
    b1 <- D/(dr^2) * (GLUCOSE[2] - GLUCOSE[1]) + D/(2*r[1]*dr) * (GLUCOSE[2] - GLUCOSE[1])
    b2 <- D/(dr^2) * (GLUCOSE[nr-1] - GLUCOSE[nr]) + D/(2*r[nr]*dr) * (GLUCOSE[nr-1] - GLUCOSE[nr])
    dGLUCOSE = c(b1, dGLUCOSE, b2)

    GlucoseUptake <- uptakeRate * GLUCOSE * MICROBE 
    MicrobeGrowth <- GlucoseUptake * CUE 
    MicrobeMortality <- mortalityRate * MICROBE
    
    dGLUCOSE = dGLUCOSE - GlucoseUptake + MicrobeMortality + base_release
    dMICROBE = MicrobeGrowth - MicrobeMortality
    
    return(list(c(dGLUCOSE, dMICROBE)))
  })
}
#------------------------------------------------------------------------------------------------------------------------
# Time
dt <- 3600/4  # time step (seconds)
days <- 10    # number of days to run the simulation
duration <- days * 86400  # days to seconds
nt <- duration/dt + 1     # number of observations
time <- seq(0, duration, dt)  # seconds

# Space (1-D)
r0 <- 0.03    # where root meets soil (mm)
distance <- 30  # distance away from center (mm)

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
#--------------------------------------------------------------------------------------------------------------------------
# Model parameters and initial conditions
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
           SOC_mineralization = SOC_mineralization)

base_release <- 0.1

# Initial conditions
Glucose <- rep(0, N)
Microbes <- rep(0.1, N)  # Adjust initial microbe concentration if needed

state <- c(Glucose, Microbes)

# Solve the model
model <- RImodel(0, state, parms, distance, N)
result <- ode.1D(y = state, times = time, func = RImodel, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time <- result[, 1]
GLUCOSE_result <- result[, 2:(N + 1)]
MICROBE_result <- result[, (N + 2):(2*N + 1)]

# Print the simulation results for glucose and microbe concentrations
print(result)

#################################################################
# Time series plot for Glucose and Microbe populations
plot(time, GLUCOSE_result[, 1], 
     type = "l", 
     xlab = "Time (s)", 
     ylab = "Glucose Concentration")
lines(time, MICROBE_result[, 1], col = "red")
legend("topright", legend = c("Glucose", "Microbes"), col = c("black", "red"), lty = 1)

# Equilibrium Microbe Plot
plot(time, cum_microbes_mass_change, 
     type = "l",
     xlab = "Time (s)", 
     ylab = "Cumulative Microbes Change")
abline(h = 0, col = "red", lty = 2)

# Total Mass Comparison Plot
plot(time, glucose_total_mass, 
     type = "l", 
     xlab = "Time (s)", 
     ylab = "Total Glucose Mass")
lines(time, microbes_total_mass, col = "red")
legend("topright", legend = c("Glucose", "Microbes"), col = c("black", "red"), lty = 1)

# SOC Mineralization Comparison Plot
plot(time, SOC_mineralization_total, 
     type = "l", 
     xlab = "Time (s)", 
     ylab = "Total SOC Mineralization")
lines(time, SOC_mineralization_zero, col = "red")
legend("topright", legend = c("Given SOC", "SOC = 0"), col = c("black", "red"), lty = 1)

# Microbe Equilibrium Comparison Plot
plot(time, MICROBE_result[, 1], 
     type = "l", 
     xlab = "Time (s)", 
     ylab = "Microbe Population")
lines(time_modified, MICROBE_result_modified[, 1], col = "red")
legend("topright", legend = c("Original Run", "Modified Initial Conditions"), col = c("black", "red"), lty = 1)
