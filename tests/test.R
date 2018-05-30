library(ggplot2)
library(FARg)
library(FARallnat)


model <- "cnrm"
# load data from the package
data(list=model)
# formating data, e.g passing from temperature to anomalie, keep only hist
# and rcp runs
mdata <- format_data(get(model))
# keep only run that posses rcp simulations
if(model != "obs") mdata <- select_continuous_run(mdata)
xp <- 1.6
xvar <- "eur_tas"
yvar <-  "eur_tas"
tvar <- "year"

compute_and_plot_far <- function(mdata, y="y", x="x", time="time", xp=1.6 , R=3, stat_model=gauss_fit, ci_p=0.9, ...){
  ans <- compute_far_simple.default(mdata,
                                    y = yvar, x = xvar, time = tvar,
                                    xp = xp, R  = R, ci_p = ci_p,
                                    stat_model = stat_model, ...) 
  ellipsis <- list(...)
  ellipsis <- sapply(ellipsis, as.character)
  if(length(ellipsis) > 0 ){
    ellipsis <- paste(names(ellipsis), ellipsis, sep = "_", collapse = "_")
  } 
  else{
    ellipsis <- ""  
  }
  pdf_name <- paste("FAR", "stat_model", deparse(substitute(stat_model)), "y", y, "x", x, "xp", xp, "R", R,  ellipsis, sep = "_")
  pdf(file = paste(pdf_name, ".pdf", sep = ""))
  
  oridat <- data.frame(year = mdata[, tvar],
                       x = mdata[, xvar],
                       y = mdata[, yvar])
  print(head(oridat))
  merged_res <- ans$ic_far
  p_dx <- default_plot(subset(merged_res, param %in% c("x_all", "x_ant", "x_nat"))) +
    geom_hline(aes_q(yintercept=xp), linetype=2) + 
    geom_point(data = oridat, aes(x=year, y=x), color="black",size=0.9, alpha=0.7) 
  p_dx <- reorder_layer(p_dx)
  plot(p_dx)
  
  p_mu <- default_plot(subset(merged_res, param %in% c("mu_all", "mu_ant", "mu_nat", "threshold_all", "threshold_ant", "threshold_nat"))) +
    geom_hline(aes_q(yintercept=xp), linetype=2) + 
    geom_point(data = oridat, aes(x=year, y=y), color="black",size=0.9, alpha=1) 
  p_mu <- reorder_layer(p_mu)
  plot(p_mu)
  
  p_p <- default_plot(subset(merged_res, param %in% c("p_all", "p_ant", "p_nat"))) + facet_grid(param ~ ., scales = "free")+
    coord_trans(y="sqrt")
  plot(p_p)
  
  plot_far(merged_res, "al")
  dev.off()
  return(ans)
}
compute_and_plot_far(mdata, y = yvar, x = xvar , time = tvar, xp=1.6 , R=3, stat_model=gauss_fit, ci_p=0.9)
compute_and_plot_far(mdata, y = yvar, x = xvar , time = tvar, xp=1.6 , R=3, stat_model=gpd_fit, ci_p=0.9, qthreshold = .9)
