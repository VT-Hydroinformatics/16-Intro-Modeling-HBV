# HBV (Hydrologiska byrÂns vattenavdelning) model 
# Code for initializing HBV hydrologic simulation model
#
# References:
# Bergstrˆm, S., & Forsman, A. (1973). Development of a conceptual 
#     deterministic rainfall-runoff model. Nordic Hydrology, 4(3), 147-170.
# Bergstrˆm, S. (1995) The HBV model (Chapter 13, pp. 443-476), in: 
#     Singh, V.P. (ed.) Computer models of watershed hydrology, Water 
#     Resources Publications, Highlands Ranch, Colorado, U.S.A., 1130 pp.
# Seibert, J. (2002) HBV light version 2, User's manual, 
#     Environmental Assessment, SLU, Uppsala
#
# Matlab code written by K. McGuire, March, 2012. 
# Thanks to J. Seibert for sharing parts of HBV Light
# Code written by K. McGuire, March 2012
# Translated to R by J.P. Gannon, March 2021

# HBV (Hydrologiska byrÂns vattenavdelning) model 
#   (version does not include elevation/vegetation zones)
#
# INPUTS:
#   pars    = parameters (listed below in order)
#   P       = precipitation daily time series, mm/d
#   T       = mean daily temperature, deg C
#   PET     = estimate of daily potential evapotranspiration, mm/d
#   routing = when this term = 1, then triangular routing is invoked, or
#             routing = 0 
#
# OUTPUTS:
#   q       = simulated streamflow, mm/d
#   qs      = simulated shallow flow (e.g., SOF, macropores, etc.), mm/d
#   qi      = simulated interflow (i.e., from upper GW storage), mm/d
#   qb      = simulated baseflow (i.e., from lower GW storage), mm/d
#   Storage = simulated total storage (= soil + GW1 + GW2), mm
#   SWE     = simulated snow water equivalent, mm
#   AET     = simulated actual evapotranspiration, mm/d
#   SF      = simulated snowfall (includes evap and catch efficiency), mm/d
#   S1      = simulated storage from the upper GW reservoir, mm
#   S2      = simulated storage from the lower GW reservoir, mm
#   soil    = simulated storage in soil, mm
#   w       = simulated water inputs to soil (= rainfall + snowmelt), mm/d


HBV <- function(pars,P,Temp,PET,routing){
  # [q,qs,qi,qb,Storage,SWE,AET,SF,S1,S2,soil,w] = HBV(pars,P,T,PET)

  ## PARAMETERS
  FC    <- pars[1]  #Max soil moisture storage, field capacity
  beta  <- pars[2]  #Shape coefficient governing fate of water input to soil moisture storag  
  LP    <- pars[3]  #Threshold for reduction of evap
  SFCF  <- pars[4]  #Snowfall correction factor
  TT    <- pars[5]  #Threshold temperature
  CFMAX <- pars[6]  #Degree-day factor
  CFR   <- 0.05     #Usually FIXED  Refreezing coefficient (use default value of 0.05)
  CWH   <- 0.1      #Usually FIXED  Water holding capacity (use default value of 0.1)
  k0    <- pars[7]  #Recession constant (upper storage, near surface)
  k1    <- pars[8]  #Recession constant (upper storage)
  k2    <- pars[9]  #Recession constant (lower storage)
  UZL   <- pars[10] #Threshold for shallow storage
  PERC  <- pars[11] #Percolation, max flow from upper to lower storage
  MAXBAS<- pars[12] 
  
  ## INITIALISE VARIABLES
  AET     <- rep(0, length(P))   #Actual evap, mm/d
  SF      <- rep(0, length(P))   #Snowfall, mm
  R       <- rep(0, length(P))   #Recharge, mm/d
  soil    <- rep(0, length(P))   #Soil storage, mm
  SWE     <- rep(0, length(P))   #SWE, snow water equivalent, mm
  w       <- rep(0, length(P))   #water input = snowmelt + rainfall, mm/d
  Storage <- rep(0, length(P))   #Total storage, mm
  S1      <- rep(0, length(P))   #Upper zone storage, mm
  S2      <- rep(0, length(P))   #Lower zone storage, mm
  Q_STZ   <- rep(0, length(P))   #unrouted shallow flow = Saturated overland flow or other rapid process, mm/d
  Q_SUZ   <- rep(0, length(P))   #unrouted upper zone flow = interflow, mm/d
  Q_SLZ   <- rep(0, length(P))   #unrouted lower zone flow = groundwater flow, mm/d
  Qgen    <- rep(0, length(P))   #Streamflow sources not routed through channel network
  
  SWE[1] <- 0 
  SUZ    <- 0       #Initial upper zone storage
  SLZ    <- 0       #Initial lower zone storage
  SP     <- 152.4   #initial value for simulated snowpack 
  WC     <- 0       #Initial liquid water in snowpack 
  SM     <- FC      #Initial soil storage content
  
  
  ## --------- TIME LOOP ------------
  for (t in 2:length(P)){
    
    ## SNOW
    inc<- 0 
    
    if (SP > 0){  
      if (P[t] > 0){  
        if (Temp[t] > TT){
          WC <- WC + P[t] 
        }else{
          SP <- SP + P[t] * SFCF 
        }
      }
      
      if (Temp[t] > TT){
        melt <- CFMAX * (Temp[t] - TT)   #Equation for melting rate of snowpack, but not more than snow is available
        if (melt > SP){ 
          inc<- SP + WC 
          WC <- 0 
          SP <- 0 
        }else{
          SP <- SP - melt        #Total simulated snow is the simulated snowpack for the day before minus the melted snow
          WC <- WC + melt        #Total water amount is increased by melted snow
          if (WC >= CWH * SP){   #The snowpack can retain as much as 10# of its water equivalent, but not more
            inc<- WC - CWH * SP  #If there is more liquid water, this goes to runoff (note:if there is no snowpack all water will go to catchment input)
            WC <- CWH * SP 
          }
        }
      }else{
        refreeze <- CFR * CFMAX * (TT - Temp[t])   #Refreezing of meltwater that occurs if T < TT, but not more than liquid water is available
        if (refreeze > WC){
          refreeze <- WC  
        }
        SP <- SP + refreeze   #The total simulated snowpack includes refreeze of that day.
        WC <- WC - refreeze   #The total simulated amount of water in the snowpack is the water in the snowpack the day before minus the water refrozen the same day
        SF[t] <- P[t] * SFCF  #Snowfall
      }
    }else{
      if (Temp[t] > TT){
        inc<- P[t]   #If too warm, input is rain
      }else{
        SP <- P[t] * SFCF 
      }
    }
    SWE[t] <- SP  + WC      #SWE, snow water equivalent, mm
    
    
    ## SOIL MOISTURE 
    rQ <- 0 
    old_SM <- SM 
    if (inc> 0){
      if (inc< 1){
        y = inc
      }else{
        m = floor(inc)  #loop through 1 mm increments
        y = inc- m 
        for (i in 1:m){   #Loop for adding input to soil 1 mm at a time to avoid instability
          dQdP = (SM/FC)^beta     #Partitioning function btw soil moisture storage and recharge
          if (dQdP > 1)  dQdP = 1
          SM <- SM + 1 - dQdP 
          rQ <- rQ + dQdP 
        }
      }
      dQdP <-(SM/FC)^beta 
      if (dQdP > 1)  dQdP = 1
      SM = SM + (1 - dQdP)*y     #Soil moisture update
      rQ = rQ + dQdP*y           #Recharge
    }
    
    mean_SM <- (SM + old_SM)/2     #average previous and current soil moisture for estimating AET
    if (mean_SM < LP*FC){
      AET[t] <- PET[t]*mean_SM/(LP*FC)  #AET estimated as simple linear function with storage deficits
    }else{
      AET[t] <- PET[t] 
    }
    
    if (SP + WC > 0){
      AET[t] <- 0       #No evap if snow present (SFCF accounts for this and catch error)
    }
    SM = SM - AET[t]      #Update soil moisture with AET flux
    
    if (SM < 0){
      soil[t] <- 0  
      SM <- 0 
    }
    soil[t] <- SM 
    R[t] <- rQ 
    w[t] <- inc 
    
    
    ## GW STORAGE
    SUZ <- SUZ + rQ        #add recharge to upper storage
    if (SUZ - PERC < 0){   #if percolation > upper storage, all water moves ot lower storage
      SLZ <- SLZ + SUZ 
      SUZ <- 0 
    }else{
      SLZ <- SLZ + PERC    #increase lower storage with percolation
      SUZ <- SUZ - PERC 
    }
    
    if (SUZ < UZL){
      Q_STZ[t] <- 0      #if upper storage < UZL, then there is no shallow flow
    }else{
      Q_STZ[t] = (SUZ-UZL)*k0     #shallow flow calc
    }
    
    Q_SUZ[t] <- SUZ*k1                        #flow from upper storage
    Q_SLZ[t] <- SLZ*k2                        #flow from lower storage
    SUZ      <- SUZ - Q_SUZ[t] - Q_STZ[t]     #Upper storage state update
    SLZ      <- SLZ - Q_SLZ[t]                #Lower storage state update
    S1[t]    <- SUZ 
    S2[t]    <- SLZ 
    
    ##DISCHARGE COMPONENTS
    Qgen[t] <- Q_STZ[t] + Q_SUZ[t] + Q_SLZ[t] #Streamflow sources not routed through channel network
    
    ## OUTPUT STORAGE
    Storage[t] <- S1[t] + S2[t] + soil[t] #Total storage
  }
  ##--------- END TIME LOOP ------------
  
  ## ROUTING (for small catchments...this component is less important)
  
  if (routing == 1){ 
    step <- 0.005  #du or integration step
    
    i <- (0:step:MAXBAS) 
    h <- rep(0, length(i)) 
    
    j <- which(i < MAXBAS/2) 
    h(j) <- step * (i[j] * 4 / MAXBAS ^ 2)  #constructs the triangular weighting function
    
    j <- which(i >= MAXBAS/2) 
    h(j) = step * (4 / MAXBAS - i[j] * 4 / MAXBAS ^ 2)   #constructs the triangular weighting function
    
    #allow base of weighting function to be noninteger, adjusts for extra weights for the last day
    if (MAXBAS %% 1 > 0){          
      I <- (1:length(i) / MAXBAS-1:length(i)) 
      I <- I:length(i)
    }else{
      I <- (1:length(i) / floor(MAXBAS)-1:length(i)) 
    }
    
    MAXBAS_w <- rep(0, length(I)) 
    
    #integration of function
    for (k in 2:length(I)){  
      MAXBAS_w[k] = sum(h[floor(I[k-1]):floor(I[k])]) 
    }
    #make sure integration sums to unity for mass balance
    MAXBAS_w <- MAXBAS_w[2:end] / sum(MAXBAS_w[2:end])  
    
    # ROUTING OF DISCHARGE COMPONENTS
    qs <- conv(Q_STZ,MAXBAS_w)  
    qs <- qs[1:length(P)] #routed shallow flow = Saturated overland flow or other rapid process
    qi <- conv(Q_SUZ,MAXBAS_w)  
    qi <- qi[1:length(P)] #routed upper zone flow = interflow
    qb <- conv(Q_SLZ,MAXBAS_w)  
    qb <- qb[1:length(P)] #routed lower zone flow = baseflow
    q  <- conv(Qgen,MAXBAS_w)  
    q  <- q[1:length(P)]  #total flow routed to catchment outlet
  }else{
    
    ## NO ROUTING
    qs <- Q_STZ   #not routed shallow flow = Saturated overland flow or other rapid process
    qi <- Q_SUZ   #not routed upper zone flow = interflow
    qb <- Q_SLZ   #not routed lower zone flow = baseflow
    q  <- Qgen    #total flow routed to catchment outlet
  }

  ##Create tibble for output
     tibble(q,       
     qs      ,
     qi      ,
     qb      ,
     Storage ,
     SWE     ,
     AET     ,
     SF      ,
     S1      ,
     S2      ,
     soil    ,
     w       )

}
