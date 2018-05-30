#######################################################################################
# EBM
#######################################################################################
# Based on Geoffroy et al. (2013)

# Two-box model for the climate system (ocean-atmosphere).
# c  dT/dt   = F - lamb T - gamm (T - T_0)
# c0 dT_0/dt = gamm (T - T_0)
# Based on Held et al. (2010)

#-- Arguments :
# forcing : step - linear - stabilization
# F_infty : amplitude of the forcing
# c       : atm/upper ocean heat capacity
# c0      : deep ocean heat capacity
# lamb    : feedback parameter
# gamm    : ocean heat uptake coefficient
# N       : number of points simulated   

held_model  <- function (FF, model="CNRM"){
  data(ebm_params,  envir = environment())
  model_param <- subset(ebm_params, models_name == model)
  list2env(model_param, environment())
  #-0- Forcing functions
  N <- length(FF)
  dt <- 1; #-- timestep (year)

  #-1- Numerical solutions (explicit scheme)
  T <- numeric(N+1);
  To <- numeric(N+1);
  H <- numeric(N+1);
  T[1] <- 0;
  To[1] <- 0;
  for(n in 2:(N+1)){
    T[n]  <-  (T[n-1] + dt/c*(FF[n-1] - lamb*T[n-1] - gamm*(T[n-1]-To[n-1])));
    To[n] <-  (To[n-1] + dt/c0*(gamm*(T[n-1]-To[n-1])));
    H[n]  <-  gamm*(T[n] - To[n]);
  }
  ans <- cbind(T, To, H)
  colnames(ans) <- c("T", "To", "H")
  ans
}


