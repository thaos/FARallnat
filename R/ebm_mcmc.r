##########################################################################
# bayesian : likelihood, prior and posterior.

make_likelihood <- function(y, y_year){
  data(FF)
  FF <- FF$FF
  f_ghg <- FF$ghg
  f_nat <- FF$nat
  f_others <- FF$others
  i_ebm <- fmatch(y_year, FF$year)
  likelihood <- function(param){
    int = param[1]
    sf_ghg = param[2]
    sf_nat = param[3]
    sf_others = param[4]
    c = param[5]
    c0 = param[6]
    gamm = param[7]
    sd = param[8]
    if(sd < 0) return(Inf)
    pred <- hmodel_sf(f_ghg, f_nat, f_others, int, sf_ghg, sf_nat, sf_others,  c, c0, lamb=1, gamm)[i_ebm]
    singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)
    sumll = sum(singlelikelihoods)
    return(sumll)   
    }
}
# likelihood <- make_likelihood(y=tas_cnrm$gbl_tas, y_year=tas_cnrm$year)

# param_cnrm <- c("int" = mean(head(tas_cnrm$gbl_tas, n=30*5)),
#                 "sf_ghg" = 1,
#                 "sf_nat" =1,
#                 "sf_others" = 1,
#                 "c"=parameters[3, 5],
#                 "c0"=parameters[3, 6],
#                 "gamm"= parameters[3, 4],
#                 "sd"=sd(head(tas_cnrm$gbl_tas, n=30*5)))
# likelihood(param_cnrm)
# 

# Prior distribution
make_prior <- function(y){
  data(FF)
  sf_sig <- FF$sig
  sig_ghg <- sf_sig["ghg"]
  sig_nat <- sf_sig["nat"]
  sig_others <- sf_sig["others"]
  mean30 <- mean(head(y, n=30))
  sd30 <- sd(head(y, n=30))
  p_avg <- apply(parameters[, -1], 2, mean)
  p_sd <- apply(parameters[, -1], 2, sd)
  c_mean <- p_avg["c"]
  c_sd <- p_sd["c"]
  c0_mean <- p_avg["c0"]
  c0_sd <- p_sd["c0"]
  gamm_mean <- p_avg["gamm"]
  gamm_sd <- p_sd["gamm"]
  start_value <- c(mean30, 1, 1, 1, c_mean, c0_mean, gamm_mean, sd30)
  prior <- function(param){
    int = param[1]
    sf_ghg = param[2]
    sf_nat = param[3]
    sf_others = param[4]
    c = param[5]
    c0 = param[6]
    gamm = param[7]
    sd = param[8]
    int_prior = dnorm(int, mean=mean30,  sd = sd30, log = T)
    sf_ghg_prior = dnorm(sf_ghg, mean=1,  sd =sig_ghg , log = T)
    sf_nat_prior = dnorm(sf_nat, mean=1,  sd =sig_nat , log = T)
    sf_others_prior = dnorm(sf_others, mean=1,  sd =sig_others , log = T)
    c_prior = dnorm(c, mean=c_mean,  sd=c_sd , log = T)
    c0_prior = dnorm(c0, mean=c0_mean,  sd=c0_sd , log = T)
    gamm_prior = dnorm(gamm, mean=gamm_mean,  sd=gamm_sd , log = T)
    sd_prior = dunif(sd, min=0, max=sd30*3, log = T)
    return(int_prior +
           sf_ghg_prior +
           sf_nat_prior +
           sf_others_prior +
           c_prior +
           c0_prior +
           gamm_prior +
           sd_prior)
  }
  attr(prior, "start_value") <- start_value
  prior
}
# prior <- make_prior(y=tas_cnrm$gbl_tas)


posterior <- function(param){
     return (likelihood(param) + prior(param))
}

######## Metropolis algorithm ################
 
proposalfunction <- function(param){
      ans <- rnorm(8,mean = param, sd= c(0.05, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01))
      if(ans[8] < 0) ans[8] <- abs(ans[8])
      ans
}
 
run_metropolis_MCMC <- function(startvalue, iterations){
  chain = array(dim = c(iterations+1,length(startvalue)))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    cat(i, "proposal : ", proposal[c(2, 8)], "\n")
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

continue_metropolis_MCMC <- function(chain, iterations){
  chain <- rbind(chain, run_metropolis_MCMC(tail(chain, 1), iterations))
}

 
# startvalue = attr(prior, "start_value")
# print(system.time({chain = run_metropolis_MCMC(startvalue, 100000)}))
# print(system.time({chain = continue_metropolis_MCMC(chain, 100000)}))
# print(system.time({chain = run_metropolis_MCMC(startvalue, 1000000)}))
 
# burnIn = 70000
# acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

### Summary: #######################

get_param <- function (chain, burnIn=1){
  ans <- apply(chain, 2, median)
  ans
}

predict.mcmc <- function(param, f_ghg, f_nat, f_others){
    int = param[1]
    sf_ghg = param[2]
    sf_nat = param[3]
    sf_others = param[4]
    c = param[5]
    c0 = param[6]
    gamm = param[7]
    # gamm = 0.48
    sd = param[8]
    pred <- hmodel_sf(f_ghg, f_nat, f_others, int, sf_ghg, sf_nat, sf_others,  c, c0, lamb=1, gamm)
    pred
}
# param <- get_param(chain, burnIn=150000)

plot_mcmc <- function(chain, burnIn=1){
  print(acceptance  <-  1-mean(duplicated(chain[-(1:burnIn),])))
  chain <- chain[-(1:burnIn), ]
  estim <- get_param(chain, burnIn=burnIn)
  estim <- data.frame(estim=estim)
  estim$variable <- c("int", "sf_ghg", "sf_nat", "sf_others", "c", "c0", "gamm", "sd") 
  chain <- as.data.frame(chain)
  names(chain)= estim$variable
  chain$id <- 1:nrow(chain)
  chain_m <- melt(chain, id=c("id"))
  p_chain <- ggplot(data=chain_m, aes(x=id, y=value))+
  geom_hline(data=estim, mapping=aes(yintercept=estim), color="blue") +
  geom_point()+
  facet_wrap(~variable, scales="free")
  p_hist <- ggplot(data=chain_m)+
  geom_histogram(aes(x=value))+
  geom_vline(data=estim, mapping=aes(xintercept=estim), color="blue") +
  facet_wrap(~variable, scales="free")
  grid.arrange(p_hist, p_chain, ncol=1)
}
# plot_mcmc(chain)
# plot_mcmc(chain, burnIn=100000)

