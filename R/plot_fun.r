#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @import scales
#' @importFrom grDevices dev.off hcl pdf
#' @importFrom graphics plot
#' @importFrom stats aggregate dnorm dunif lm median optim predict quantile rnorm runif sd time
#' @importFrom utils data head tail


# Plot results of computer far.
# From the CI data frame plot the time evoltuion of the confidence intervals for
# a given parameter
# Creates one panel for each value of the column method
plot_pannel_boot_time <- function(ic_df, param="FAR", main="FAR(t1)"){ 
  env <- environment()
  ic_subset <- ic_df[ic_df$param == param, ]
  ggplot(ic_subset,  aes_string(x="time"))+
  ggtitle(main)+
  geom_point(aes_string(x="time",  y="Estim",  group="method",  color="method"), size=1)+
  geom_line(aes_string(x="time",  y="Estim",  group="method",  color="method"), size=1)+ 
  geom_ribbon(aes_string(x="time",  ymin="IC_inf",  ymax="IC_sup",  group="method",  fill="method"),  alpha=0.2)+ 
  facet_wrap(~method)# +
  # coord_cartesian( xlim = xlim, ylim = ylim)+
  # theme(legend.position = "bottom")
}

# a given parameter
# Make a plot to visualize the evolution of one one parameter 
#' \code{plot_far} Plot the times series of one model parameters and its
#' confidence intervales
#' @param ic_df a data_frame contaning the parameter estimates and their confidence
#' intervals
#' @param param the name of the parameter to be plotted 
#' @param main the title of the plot
#' @param alpha the transparency of the shading for the confidence intervals
#' @examples
#' library(FARg)
#' model <- "cnrm"
#' #load data from the package
#' data(list=model)
#' # formating data, e.g passing from temperature to anomalie, keep only hist
#' # and rcp runs
#' mdata <- format_data(get(model))
# keep only run that posses rcp simulations
#' if(model != "obs") mdata <- select_continuous_run(mdata)
#' ans <- compute_far_simple.default(mdata,
#'  y="eur_tas", x="gbl_tas", time="year",
#'   xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' plot_boot_time(ans$ic_far, param="p_all", main="p_all") 
#' @export
plot_boot_time <- function(ic_df, param="FAR", main="", alpha=0.2){ 
  env <- environment()
  ic_subset <- ic_df[ic_df$param == param, ] 
  ggplot(ic_subset,  aes_string(x="time"))+
    ggtitle(main)+
    geom_point(aes_string(x="time",  y="Estim",  group="method",  color="method"), size=1)+
    geom_line(aes_string(x="time",  y="Estim",  group="method",  color="method"), size=1)+ 
    geom_ribbon(aes_string(x="time",  ymin="IC_inf",  ymax="IC_sup",  group="method",  fill="method"),  alpha=alpha)+ 
  # coord_cartesian( xlim = xlim, ylim = ylim)+
  # theme(legend.position = "bottom")
  theme(legend.position = "right")
}

#' generate a vactor of n equally spaced (in terms of hue)  colors.
#' 
#' \code{gg_color_hue} generate a vector of n equally spaced (in terms of hue)
#' colorsi for ggplot2 plots.
#' @param n the number of color wanted
#' @return a vector of colors
#' @examples
#' gg_color_hue(10)
#' @export
# generates a set of n colors for ggplot2 plots
gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# Make a double y-axis plot with RR and the FF
#' \code{plot_far} create a symetric double y-axis plot for the evolution FAR and the RR
#' @param ic_mat a data_frame contaning the RR estimates and their confidence
#' intervales and the rescaled value of RR (either RRal or RRel, see axis_trans)
#' @param axis_trans the choice of doublei : axis = "al" for a symetric
#' atan(log(y)) axis ; axis_trans = "el" for a symetrized FAR axis 
#' @param main the title of the plot
#' @param xlim a vector of two values for the limits of the x-axis
#' @param ylim a vector of two values for the limits of the y-axis
#' @param col a vector of colors for custom colors, should be the same length
#' as the number of different values of the ic_mat method column.
#' @examples
#' library(FARg)
#' model <- "cnrm"
#' #load data from the package
#' data(list=model)
#' # formating data, e.g passing from temperature to anomalie, keep only hist
#' # and rcp runs
#' mdata <- format_data(get(model))
# keep only run that posses rcp simulations
#' if(model != "obs") mdata <- select_continuous_run(mdata)
#' ans <- compute_far_simple.default(mdata,
#'  y="eur_tas", x="gbl_tas", time="year",
#'   xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' plot_far(ans$ic_far, axis_trans="al", col=gg_color_hue(2))
#' @export
plot_far <- function(ic_mat, axis_trans, main="", xlim, ylim, col=NULL){
  ticks=c(0, 1/100, 1/10, 1/5, 1/3, 1/2, 1, 2, 3, 5, 10, 100, Inf)
  if(missing(xlim)) xlim <- range(ic_mat$time)
  if(missing(ylim)) ylim <- c(-1.1, 1.1)
  inv <- get(paste(axis_trans, "_inv", sep=""))
  trans <- get(paste(axis_trans, "_trans", sep=""))
  breaks <- trans(ticks) 
  param <- paste(axis_trans,"RR", sep="")
  p1 <- plot_boot_time(ic_mat, param=param, main=main, alpha=0.2) +
  scale_y_continuous(name = "Relative Risk",
                     breaks = breaks,
                     labels = scales::trans_format(inv, function(x) format(x, digits=2)))+
    coord_cartesian(xlim=xlim, ylim=ylim) +
    theme(legend.position = "bottom")    
  if (!is.null(col)){ 
    p1 <- p1 + scale_color_manual(values=col)
    p1 <- p1 + scale_fill_manual(values=col)
}
  p2 <- p1 + scale_y_continuous(name="FAR",
                                breaks = breaks,
                                labels = scales::trans_format(inv, function(x) format(rrtofar(x), digits=3))) +
    coord_cartesian(xlim=xlim, ylim=ylim) +
    theme(legend.position = "bottom")    
  ggplot_dual_axis(p1,p2)
}

# plot the relative risk by pannel : each pannel correspond to a different value
# of the method column. only a sigle y-axis is available
plot_pannel_far <- function(ic_mat, axis_trans, main="", xlim, ylim, col=NULL){
  ticks=c(0, 1/100, 1/10, 1/5, 1/3, 1/2, 1, 2, 3, 5, 10, 100, Inf)
  if(missing(xlim)) xlim <- range(ic_mat$time)
  if(missing(ylim)) ylim <- c(-1.1, 1.1)
  inv <- get(paste(axis_trans, "_inv", sep=""))
  trans <- get(paste(axis_trans, "_trans", sep=""))
  breaks <- trans(ticks) 
  param <- paste(axis_trans,"RR", sep="")
  p1 <- plot_pannel_boot_time(ic_mat, param=param, main=main) +
  scale_y_continuous(name = "Relative Risk",
                     breaks = breaks,
                     labels = scales::trans_format(inv, function(x) format(x, digits=2)))+
  coord_cartesian(xlim=xlim, ylim=ylim)
  if (!is.null(col)) 
    p1 <- p1 + scale_color_manual(values=col)+ scale_fill_manual(values=col)
  p1
}

al_trans <- function(x) atan(log(x))/(pi/2)
al_inv <- function(x) exp(tan((pi/2)*x))
el_trans <- function(x) ifelse(x <= 1, x-1, 1-1/x)
el_inv <- function(x) ifelse(x > 0, 1/(1-x), x + 1)
rrtofar <- function(x) 1 - 1/x

# Confidence interval plot for results from the compute_far functions
#' \code{default_plot} Confidence interval plot for results from the compute_far functions
#' @param ic_far a data_frame contaning estimates and their confidence
#' intervales produced by the compute_far functions, i.e. the ic_far element of the returned list of the compute_far functions
#' @examples
#' library(FARg)
#' model <- "cnrm"
#' #load data from the package
#' data(list=model)
#' # formating data, e.g passing from temperature to anomalie, keep only hist
#' # and rcp runs
#' mdata <- format_data(get(model))
# keep only run that posses rcp simulations
#' if(model != "obs") mdata <- select_continuous_run(mdata)
#' ans <- compute_far_simple.default(mdata,
#'  y="eur_tas", x="gbl_tas", time="year",
#'   xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' ic_x <- subset(ans$ic_far, param %in% c("x_all", "x_ant", "x_nat"))
#' default_plot(ic_x)
#' @export
default_plot <- function(ic_far){
  ggplot(data=ic_far, aes(x=time))+
    # ggtitle(expression(paste("Covariate Estimation : ", x[t], " = ", x["t, ant"], " + ",  x["t, nat"]))) +
    geom_line(aes_string(x="time", y="Estim", group="param", color="param"), size=1) +
    geom_ribbon(aes_string(x="time", ymin="IC_inf", ymax="IC_sup", group="param", fill="param"), alpha=0.4) +
    # facet_grid(method~param) +
    facet_wrap(~method) +
    coord_cartesian(xlim=c(1850,2100))+
    # scale_color_manual(values=col)+
    # scale_fill_manual(values=col)+
    theme(legend.position = "bottom")
}

# Put the bottom layer of a ggplot on top
#' \code{reorder_layer} Put the bottom layer of a ggplot on top
#' @param gp a ggplot2 plot 
#' @examples
#' library(FARg)
#' model <- "cnrm"
#' #load data from the package
#' data(list=model)
#' # formating data, e.g passing from temperature to anomalie, keep only hist
#' # and rcp runs
#' mdata <- format_data(get(model))
# keep only run that posses rcp simulations
#' if(model != "obs") mdata <- select_continuous_run(mdata)
#' ans <- compute_far_simple.default(mdata,
#'  y="eur_tas", x="gbl_tas", time="year",
#'   xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' plot_far(ans$ic_far, axis_trans="al", col=gg_color_hue(2))
#' @export
reorder_layer  <- function(gp){
  layers <- gp$layer
  nb_layers <- length(layers)
  if(nb_layers >1) layers <- layers[c(2:nb_layers, 1)]
  gp$layers <- layers
  gp
}