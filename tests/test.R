library(devtools)
load_all("../farg")
load_all()
model <- "cnrm"
# load data from the package
data(list=model)
# formating data, e.g passing from temperature to anomalie, keep only hist
# and rcp runs
mdata <- format_data(get(model))
# keep only run that posses rcp simulations
if(model != "obs") mdata <- select_continuous_run(mdata)

ans <- compute_far_simple.default(mdata,
                                  xp = 1.6,
                                  R  = 100,
                                  y = "eur_tas",
                                  x = "gbl_tas",
                                  time = "year",
                                  ci_p=0.9) 

plot_boot_time(ans$ic_far, param = "p_all")
plot_far(ans$ic_far, "al")

default_plot <- function(subset){
  ggplot(data=subset, aes(x=time))+
    # ggtitle(expression(paste("Covariate Estimation : ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
    geom_line(aes(x=time, y=Estim, group=param, color=param), size=1) +
    geom_ribbon(aes(x=time, ymin=IC_inf, ymax=IC_sup, group=param, fill=param), alpha=0.4) + 
    # facet_grid(method~param) +
    facet_wrap(~method) +
    coord_cartesian(xlim=c(1850,2100))+
    # scale_color_manual(values=col)+
    # scale_fill_manual(values=col)+
    theme(legend.position = "bottom")
}
reorder_layer  <- function(gp){
  layers <- gp$layer
  nb_layers <- length(layers)
  if(nb_layers >1) layers <- layers[c(2:nb_layers, 1)]
  gp$layers <- layers
  gp
}
p_dx <- default_plot(subset=subset(merged_res, param %in% c("x_all", "x_ant", "x_nat"))) +
  geom_hline(aes_q(yintercept=xp), linetype=2)+
  geom_point(aes_string(y=x), color="black",size=0.9, alpha=0.1) 
p_dx <- reorder_layer(p_dx)
p_mu <- default_plot(subset=subset(merged_res, param %in% c("mu_all", "mu_ant", "mu_nat", "threshold_all", "threshold_ant", "threshold_nat"))) +
  geom_hline(aes_q(yintercept=xp), linetype=2)+
  geom_point(aes_string(y=y), color="black",size=0.9, alpha=0.1) 
p_mu <- reorder_layer(p_mu)
p_p <- default_plot(subset=subset(merged_res, param %in% c("p_all", "p_ant", "p_nat"))) +
  coord_trans(y="sqrt")
subset_q <- subset(merged_res, param %in% c("q_all", "q_ant", "q_nat"))