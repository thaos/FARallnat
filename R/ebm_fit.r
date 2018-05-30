#' @importFrom fastmatch fmatch 

# EBM model function
hmodel <- function(FF, c, c0, lamb, gamm){
  N <- length(FF)
  dt <- 1; #-- timestep (year)
  #-1- Numerical solutions (explicit scheme)
  T <- numeric(N+1);
  To <- numeric(N+1);
  # H <- numeric(N+1);
  T[1] <- 0;
  To[1] <- 0;
  for(n in 2:(N+1)){
    T[n]  <-  (T[n-1] + dt/c*(FF[n-1] - lamb*T[n-1] - gamm*(T[n-1]-To[n-1])));
    To[n] <-  (To[n-1] + dt/c0*(gamm*(T[n-1]-To[n-1])));
    # H[n]  <-  gamm*(T[n] - To[n]);
  }
  # ans <- cbind(T, To, H)
  # colnames(ans) <- c("T", "To", "H")
  # ans
  T[-1]
}

# EBM model function with scaling factors
# data(cnrm)
# data(FF)
# data(ebm_params)
# tas_cnrm <- merge(cnrm, FF$FF, by="year")
# with(tas_cnrm, plot(hmodel_sf(ghg, nat, others, 288, 1, 1, 1, ebm_params[12, "c"], ebm_params[12, "c0"], ebm_params[12, "lamb"], ebm_params[12, "gamm"])))
hmodel_sf <- function(f_ghg, f_nat, f_others, int, sf_ghg, sf_nat, sf_others,  c, c0, lamb, gamm){
  f_ghg <- f_ghg * sf_ghg
  f_nat <- f_nat * sf_nat
  f_others <- f_others * sf_others
  r_ghg <- hmodel(f_ghg, c, c0, lamb, gamm)
  r_nat <- hmodel(f_nat, c, c0, lamb, gamm)
  r_others <- hmodel(f_others, c, c0, lamb, gamm)
  r_all <- r_ghg + r_nat + r_others + int
}


# EBM function with scaling factor used for the ordinary least squares optimization 
hmodel_sf_optim <- function(par, obs_df, forcings_df){
  # print(par)
  hsim <- hmodel_sf(forcings_df$ghg, forcings_df$nat, forcings_df$others, int= par[1], sf_ghg=par[2], sf_nat=par[3], sf_others=par[4],  c=par[5], c0=par[6], lamb=1, gamm=par[7])
  hsim = hsim[fastmatch::fmatch(obs_df$year, forcings_df$year)]
  # with(obs_df, plot(year, obs))
  # with(obs_df, lines(year, hsim, col="red"))
  ans <- sum((obs_df$obs - hsim)^2)
  # print(ans)
  ans
}

# Fit the EBM to the observations by ordinary least squares
# data(cnrm)
# ebm_fit <- fit_ebm(cnrm, time="year", obs="gbl_tas")
fit_ebm <- function(obs_df, time="year", obs="gbl_tas"){
  data(FF,  envir = environment())
  data(ebm_params,  envir = environment())
  forcings_df <- FF$FF
  obs_df <- obs_df[, c(time, obs)]
  names(obs_df)  <- c("year", "obs")
  init=c(mean(head(obs_df$obs, n=30)), 1, 1, 1, unlist(ebm_params[12, c("c", "c0", "gamm")]))
  names(init) <- c("int", "sf_ghg", "sf_nat", "sf_others", "c", "c0", "gamm")
  print(system.time({ans <- optim(par=init, fn=hmodel_sf_optim, obs_df=obs_df, forcings_df=forcings_df, control=list(reltol=0.0001))}))
  attr(ans, "class") = "ebm_fit"
  ans
}

# predict the EBM response from an EBM fit and a new set of forcings
# data(cnrm)
# data(FF)
# data(ebm_params)
# tas_cnrm <- merge(cnrm, FF$FF, by="year")
# ebm_fit <- fit_ebm(cnrm, time="year", obs="gbl_tas")
# newdata <- unique(tas_cnrm[, c("year", "ghg", "nat", "others")])
# newdata_nat <- newdata
# newdata_nat$ghg <- 0
# newdata_nat$others <- 0
# newdata_ant <- newdata
# newdata_ant$nat<- 0
# r_all <- predict.ebm_fit(ebm_fit, newdata)
# r_ant <- predict.ebm_fit(ebm_fit, newdata_ant)
# r_nat <- predict.ebm_fit(ebm_fit, newdata_nat)
predict.ebm_fit <- function(object, newdata, ...){
  par <- object$par
  ans <- hmodel_sf(newdata$ghg, newdata$nat, newdata$others, int=par[1], sf_ghg=par[2], sf_nat=par[3], sf_others=par[4],  c=par[5], c0=par[6], lamb=1, gamm=par[7])
  ans
}


# fit the EBM model on all bootstrap samples and predict the ALL, ANT and NAT
# components for each samples
boot_ebm_allnat <- function(psamples){
  data(FF,  envir = environment())
  l_ebm <- lapply(psamples$bsamples, function(data) (fit_ebm(obs_df=data, obs="x", time="time")))
  l_ebm_an <- lapply(l_ebm, predict_gno, newdata=FF$FF)
  l_ebm_an_origin <- mapply(function(bsample, ebm_an)  merge(x=bsample, y=ebm_an, by.x="time", by.y="year"), bsample=psamples$osamples, ebm_an=l_ebm_an, SIMPLIFY=FALSE)
  l_ebm_an <- mapply(function(bsample, ebm_an)  merge(x=bsample, y=ebm_an, by.x="time", by.y="year"), bsample=psamples$bsamples, ebm_an=l_ebm_an, SIMPLIFY=FALSE)
  list("l_ebm"=l_ebm, "l_x_an"=l_ebm_an, "l_x_an_origin"=l_ebm_an_origin)
}

