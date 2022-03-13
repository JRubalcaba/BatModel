############################################################### FLIGHT POWER
# Flight costs computations with package "afpt" Heerenbrink et al. (2015) <doi:10.1098/rspa.2014.0952>

require(afpt)

M <- exp(seq(-2, 7.5, length.out = 100))     # Range of body masses (g)
Awings <- exp(seq(-6, -1, length.out = 100)) # Wing surface areas (m2)
Ta = 25  # Ambient temperature (ºC)

Pchem <- Pmech <- v_matrix <- array(NA, dim=c(100,100))
for(i in 1:length(M)){
  for(j in 1:length(Awings)){
    
    myFlightCondition <- list(
      density = 101325 / (287.05 * (Ta+273)), # [kg/m3]
      gravity = 9.81, # [m/s2]
      viscosity = 18.63E-6, # [m2/s] 
      windSpeed = 0,# [m/s]
      windDir = 0 # [degrees]
    )
    
    myBat = Bird(
      massTotal = M[i]/1000, # (kg) total body mass
      wingSpan = exp(0.9425) * Awings[j]^0.4934, # (m) maximum wing span
      wingArea = Awings[j], # (m2) maximum wing area
      type = "bat"
    )
    minspeed = findMinimumPowerSpeed(myBat)
    v_matrix[i,j] = as.numeric(minspeed[2])
    Pmech[i,j] = as.numeric(minspeed[3])
    Pchem[i,j] = as.numeric(minspeed[4])
  }
  print(i)
}

############################################################### THERMOREGULATION

batmodel <- function(Tc,     # Core temperature (ºC)
                     M,      # Body mass (g)
                     Awings, # Wings area (m2)
                     kw,     # Thermal conductivity of the wings (Wm-1ºC-1)
                     shape,  # Shape factor (1: spherical, >1: ellipsoidal)
                     rf,     # Insulation layer depth (m)
                     emis,   # Emissivity of long-wave radiation (0-1)
                     Ta,     # Air temperature (ºC)
                     v,      # Wind speed (ms-1)
                     RH)      # Solar radiation (direct + scattered) (Wm-2)
{
  ## Basal Metabolic Rate ##
  BMR_O2 = exp(0.8508) * M^0.8013 #
  BMR = BMR_O2 / 172 # (W) 1W = 172 mLO2 h-1
  
  ## Ellipsoid geometry (Porter and Kearney 2009) ## 
  b = (3*(M*1e-6/(4*pi*shape)))^(1/3)
  c = b
  a = shape*b
  S2 = a^2*b^2*c^2 / (a^2*b^2 + a^2*c^2 + b^2*c^2) 
  V = 4/3*pi*a*b*c # Body volume (m3)
  eccentric = sqrt(a^2-c^2)/a
  Area = 2*pi*b^2 + 2*pi*((a*b)/eccentric)*asin(eccentric) # Total surface area (m2)
  
  ## Insulation layer ## 
  kins = 0.027  # Thermal conductivity fur (W m-1 ?C-1) - 0.027 (Porter and Kearney), 0.0394 (Hammel 1955)
  ao = a + rf   # inner radium + feather layer (m)
  bo = b + rf
  co = c + rf
  Vo = 4/3*pi*ao*bo*co  # External volume (m3)
  ecc_out = sqrt(ao^2-co^2)/ao
  Ao = 2*pi*bo^2 + 2*pi*((ao*bo)/ecc_out)*asin(ecc_out) # External surface (m2)

  ## Thermal radiation ## 
  sigma = 5.670367e-8  # Stefan-Boltzmann constant (W m-2 K-4)
  R = 4 * sigma * emis * (Ta+273)^3 # Thermal radiation coefficient (linear approx, Bakken 1975)
  
  ## Body convection ## 
  L = V^(1/3)   # Characteristic length (m) (Mitchell 1976)
  kf =  1.5207e-11*(Ta+273)^3 - 4.8574e-08*(Ta+273)^2 + 1.0184e-04*(Ta+273) - 3.9333e-04 # Thermal conductivity of air (W m-1 ºC-1)
  Pd = 101325   # Atmosferic pressure (Pa) 
  rho_da = Pd / (287.05 * (Ta+273))  # Density of dry air (kg m-3)
  Pv = exp(77.3450 + 0.0057 * (Ta+273) - 7235 / (Ta+273)) / (Ta+273)^8.2  # Water vapor pressure (Pa)
  RH = RH/100 # ransform relative humidity % to proportion 0-1
  x = RH * 0.62198 * Pv / Pd # Humidity ratio (mass_air_water / mas_dry_air)
  rho = rho_da * (1 + x) / (1 + 1.609 * x) # Density of moist air (Kg m-3)
  mu = 1.458e-6 * (Ta+273)^(3/2) / ((Ta+273) + 110.4)  # Air dynamic viscosity (Pa s)
  cp = 1000*(-3.46607e-20*Ta^6+9.121847e-17*Ta^5+1.079642e-13*Ta^4-5.714484e-10*Ta^3+5.773354e-07*Ta^2+8.976385e-06*Ta+1.00528623891845) # Specific heat capacity air (J Kg-1 ºC-1)
  
  Re = rho * v * L / mu     # Reynolds number
  Pr = mu * cp / kf         # Prandtl number
  Nu_forced = 2.0 + 0.6 * Re^(1/2) * (Pr)^(1/3) # Nusselt number forced convection
  
  Gr = rho^2 * 1/(Ta+273) * 9.8 * L^3 * (Tc - Ta) / mu^2  # Grashof number
  Nu_free = 2.0 + 0.6 * Gr^(1/4) * Pr^(1/3) # Nusselt number free convection
  
  Nu_total = (Nu_free^3 + Nu_forced^3)^(1/3)
  hc = Nu_total * kf/L # Convection coefficient body

  # Wing T profile
  d = exp(0.9425) * Awings^0.4934 / 2 # wing length (m) 
  Nu = 0.664 * Re^0.5 * Pr^0.33
  hc_wing = Nu * kf / d + R # wing convective heat transfer coefficient (Wm-2ºC-1; from Reichard et al. 2010)

  W = Awings / (2*d) # wing width (m)
  t = 2e-5 # wing thickness (m) (Cheney et al. 2015 https://royalsocietypublishing.org/doi/10.1098/rsif.2014.1286)
  Ac = W*d*t # cross-sectional surface area (m2)
  P = 2*W+2*t # wing perimeter (m)
  m = sqrt(hc_wing*P/(kw*Ac))
  
  Rwing = 1/2 * (cosh(m*d) + hc_wing/(m*d) * sinh(m*d)) / (sqrt(hc_wing*P*kw*Ac) * (sinh(m*d) + hc_wing/(m*d) * cosh(m*d))) # Wing resistance to covective cooling (ºC W-1)

  ## Metabolic Heat Production ## 
  kb = 0.5 + 6.14*b + 0.439 # Thermal conductivity body (Wm-1ºC-1) Porter and Kearney (2009)
  
  Rgeom = S2 / (2*V*kb) # Body resistance to heat conduction (ºC W-1)
  Rconv = 1 / (hc*Ao)   # Convective heat transfer resistance
  Rir = 1 / (R*Ao)           # IR emission resistance
  Rins = (bo - b) / (Ao * kins) # Resistance of insulation layer 
  
  Rbody = Rconv * Rir / (Rconv + Rir) + Rins 
  Rtot = (Rwing * Rbody) / (Rwing + Rbody) + Rgeom # Total resistance (ºC W-1)
  
  Q_gen = (Tc - Ta) / Rtot # Thermogenic heat (W)
  
  MR = Q_gen

  ## Model Output ## 
  output <- data.frame(MR, BMR, Ao=Ao, Atot=Ao + Awings, V, Ta, v, RH)
  return(output)
}

############################################################### Run Model

## PARAMETERS

efficiency = 0.23 # Pennycuick 2008
alpha = 0.57      # Proportion of by-product heat used for thermoregulation (Humphries & Careau 2011)
Tc = 37           # Core temperature (ºC)
kw = 20           # Thermal conductivity of the wings (W m-1 ºC-1)
shape = 1.1       # Body shape: ~1 (spheric), >1 (ellipsoidal)
rf = 0.001 * M^0.33  # Depth of layer of fur (m)
emis = 0.98       # Emissivity
Ta = 25           # Ambient temperature (ºC)

## Flight (mass-specific)
Pmech_mass <- Pmech
for(i in 1:100){
  Pmech_mass[,i] <- (Pmech[,i] / efficiency) / M
}

## Thermoregulation (mass-specific)
matQgen <- array(NA, dim=c(100, 100))
for(i in 1:100){
  v = v_matrix[,i]
  MR_mod <- batmodel(Tc=Tc, M, Awings=Awings[i], kw=kw, shape=shape, rf=rf, emis=emis, Ta=Ta, v=v, RH=50)
  
  Q_byproduct <- Pmech[,i] / efficiency - Pmech[,i]
  
  MR_real <- MR_mod$MR - alpha * Q_byproduct
  belowBMR <- which(MR_real < MR_mod$BMR)
  MR_real[belowBMR] <- MR_mod$BMR[belowBMR]
  
  matQgen[,i] <- MR_real / M
}

## Total mass-specific energy requirements
Met_rate <- matQgen + Pmech_mass

############################################################### Results
# Energy requierements in relation to wing surface area for a ~10 g bats
M[46]

plot((Pmech_mass[46,]) ~ log(Awings), type="l", col="darkblue", lwd=2,
     ylab="Mass-specific energy requierements", xlab="Wing surface area")
lines((matQgen[46,]) ~ log(Awings), col="darkred", lwd=2)
lines(Met_rate[46,] ~ log(Awings), lwd=2)

# Energy requierements in relation to body mass for a bat with wing surface area of ~119 m2
Awings[32]

plot((Pmech_mass[,32]) ~ log(M), type="l", col="darkblue", lwd=2, ylim=c(0,0.2),
     ylab="Mass-specific energy requierements", xlab="Wing surface area")
lines((matQgen[,32]) ~ log(M), col="darkred", lwd=2)
lines(Met_rate[,32] ~ log(M), lwd=2)

# Coordinate surface of energy requirements across different combinations of body mass and wing surface area

levels <- seq(-10, 2, length.out = 50)
contour(log(M_const), y=log(Awings), x=log(M), levels=levels)




