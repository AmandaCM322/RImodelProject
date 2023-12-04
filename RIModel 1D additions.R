#Load in libraries ----
library(tidyverse); library(RColorBrewer); library(deSolve)

# Function: Diffusion for Glucose and Microbes
RImodel2 <- function(time, state, parms, distance, nr0) {
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
    GlucoseUptake <- umax * MICROBE * GLUCOSE/(Kc+GLUCOSE)    
    MicrobeGrowth <- GlucoseUptake * CUE                # MicrobeGrowth represents the growth of microbes based on glucose uptake and Carbon Use Efficiency (CUE).It quantifies the increase in microbe population due to the availability of glucose, considering the efficiency of resource utilization (CUE).
    MicrobeMortality <- mortalityRate * MICROBE         # MicrobeMortality represents the rate at which microbes die off. It is calculated by multiplying the mortalityRate with the current concentration of MICROBE.
    depolymerization <-Vprime * MICROBE/(MICROBE+Km)    # Depolymerization represents the breakdown of complex molecules into simpler units by microbes.Vmaxprime is the maximum rate of depolymerization.
    # MICROBE is the microbe concentration, and km is the Michaelis constant. The rate of depolymerization increases with higher MICROBE concentration but levels off at higher concentrations due to the Michaelis constant.
    
    #Rate of change = Flux gradient + Biology 
    dGLUCOSE <-  dGLUCOSE - GlucoseUptake + f_doc * MicrobeMortality + depolymerization 
    # Microbe mortality-leading to the release of glucose back into the environment & Microbial depolymerization, breaking down complex molecules into glucose 
    dMICROBE <- MicrobeGrowth - MicrobeMortality                                      # The change in microbe population over time is influenced by: Microbial growth due to glucose uptake & Microbe mortality, representing the death of microbes 
    
    
    return(list(c(dGLUCOSE, dMICROBE))) # Returning the derivative vectors for glucose diffusion and microbe dynamics.
  })
}

# Time
dt = 7*86400 #time step (seconds)
days = 7 #number of days to run the simulation
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

#C_add <- SOC * 0.01 * pi * r[N]^2 #Add 1% of SOC in terms of glucose carbon (mg)
C_add <- 1680 * 0.05 
glu_needed <- C_add * 180.156 / 72.06 * 1000 #glucose plant production (mg)

D <- 0.5 * 10^-4 #diffusion of glucose mm2/s #sourced from Chenu and Roberson (1996) given VWC ~ 30

# Model application

#################################
# Estimate parameters for biology
#################################
secperyear = 86400*365
secperday = 86400            # relationship between eq. biomass and KM Km = Km_amp*M
basemin = 0.02/15/secperday #[ug mm-3] -+> #2percent carbon, 15years turnover 
BD = 1.1

# base mineralization 1.1 ug CO2 hr-1 kg-1 (Shelby)
basemin = 1.1/1000 * BD * 12/44/3600  # Shelby's base respiration [CO2 hr-1 kg-1]  
death = 12/secperyear  #s-1
cue  = 0.5  #carbon use efficiency | unitless

# Values to derive parameters 
#base mineralization 1.1 mg CO2 hr-1 kg-1 -> Shelby (converting into ugC mm-3 s-1)

basemin = 1.1e-3 * BD * 12/44/3600  
X = 1                        #relationship between eq. biomass and KM: Km = X*M

#Derivation of Direct parameters
doc_turn_time = 86400     # 1 day in S-1 
f_doc = 0.3               #fraction of microbial death become doc (a)
doubling_time = 30 * 86400  #Maximum microbial growth rate expressed in time for 2x
P = 0.2                   # Priming factor 2% 

c0 = basemin*doc_turn_time/ (1 - f_doc*cue)
m0 = cue*basemin/ (death * (1 - f_doc*cue)) # = aD0 

# Derived parameters and state variables equilibria
umax = log(2)/(cue * doubling_time)
Kc = umax*m0*doc_turn_time-c0

#depolymerization rate(base release under equilibrium)
Vprime = basemin / (1 - P)  

#Km calculated based on equilibrium biomass
Km = P / (1 - P) * m0 

# = 1e-6 #mm^2 s-1 
I = 0.28              #(impedance at 100% water)
Diffusion_coef = 6.3e-4 * I * 0.7  #Priesack & Kisser-Priesack (1993) 

#------------------
# Model parameters
#------------------
R <- 20     #Radius of the system
N <- 100    #Number of spatial grid points
dr <- R/N   #Spatial grid size
r <- seq(dr/2, by = dr, len = N)
ri <- seq(0, by = dr, len = N+1) #Inner radius of each grid point
dri <- dr   #Grid spacing

parms <- c(D = Diffusion_coef,
           uptakeRate = uptake, #Initial uptake rate of glucose by microbes
           CUE = cue,           #Carbon use efficiency
           umax = umax,
           f_doc,               #Fraction of microbial death become doc
           mortalityRate= death,    #Microbial mortality rate
           Km = Km,             #Michaelis-Menten constant for glucose uptake
           Kc = Kc,             #Half saturation constant for doc
           Vmaxprime = Vprime)  #Maximum depolymerization rate

# Initial conditions

#The initial glucose concentration at the first grid point is calculated as 1% of SOC (Soil Organic Carbon) concentration
Glucose = rep(c0, N) #Initialize an array of length N (number of grid points) with initial glucose concentration set to c0
Microbes = rep(m0, N) #The initial microbial concentration is set to m0

#1
Glucose [1] = c0 + C_add  #Glucose addition on & off 

#'state' is a concatenated vector containing the initial concentrations of glucose followed by microbial concentrations.
state = c(Glucose,Microbes) #Combine the arrays of initial glucose and microbial concentrations into a single state vector

# Time points for output

# Solve the model using the default banded reformulation method
model <- RImodel2(0, state, parms, distance, N)
result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

time = result[,1]

# 8 WEEKS INDIVIDUAL GLUCOSE ADDITIONS 
# Extracting results


GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]

#2
c0 = GLUCOSE_result
c0[1] = c0[1] + C_add 
m0 = MICROBE_result
state = c(c0,m0)

result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]

#3
c0 = GLUCOSE_result
c0[1] = c0[1] + C_add 
m0 = MICROBE_result
state = c(c0,m0)

result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]

#4
c0 = GLUCOSE_result
c0[1] = c0[1] + C_add 
m0 = MICROBE_result
state = c(c0,m0)

result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]

#5
c0 = GLUCOSE_result
c0[1] = c0[1] + C_add 
m0 = MICROBE_result
state = c(c0,m0)

result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")


GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]

#6
c0 = GLUCOSE_result
c0[1] = c0[1] + C_add 
m0 = MICROBE_result
state = c(c0,m0)

result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]

#7
c0 = GLUCOSE_result
c0[1] = c0[1] + C_add 
m0 = MICROBE_result
state = c(c0,m0)

result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]

#8
c0 = GLUCOSE_result
c0[1] = c0[1] + C_add 
m0 = MICROBE_result
state = c(c0,m0)

result <- ode.1D(y = state, times = time, func = RImodel2, parms = parms, distance = R, nr0 = N, dimens = N,
                 nspec = 2, names = c("GLUCOSE", "MICROBE"), method = "lsoda")

GLUCOSE_result <- result[2,2:(N + 1)]
MICROBE_result <- result[2,(N + 2):(2*N + 1)]


# Space Plots
plot(r,GLUCOSE_result)
plot(r,MICROBE_result)