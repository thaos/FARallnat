#' @import ggplot2
#' @import gridExtra
#' @import scales

# Plot results of computer far.
# From the CI data frame plot the tiem evoltuion of the confidence intervals for
# a given parameter
# Creates on panel for each value of the column method
plot_pannel_boot_time <- function(ic_df, param="FAR", main="FAR(t1)"){ 
  env <- environment()
  ic_subset <- ic_df[ic_df$param == param, ]
  ggplot(ic_subset,  aes(x=time))+
  ggtitle(main)+
  geom_point(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+
  geom_line(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+ 
  geom_ribbon(aes(x=time,  ymin=IC_inf,  ymax=IC_sup,  group=method,  fill=method),  alpha=0.2)+ 
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
#' ans <- compute_forall(compute_far.dx_ebm_fit, l_models=c("cnrm", "ipsl"), stat_model=gpd_fit, qthreshold =0.90, p=0.01, R=3)
#' plot_boot_time(ans, param="p_all", main="p_all") 
#' @export
plot_boot_time <- function(ic_df, param="FAR", main="", alpha=0.2){ 
  env <- environment()
  ic_subset <- ic_df[ic_df$param == param, ] 
  ggplot(ic_subset,  aes(x=time))+
  ggtitle(main)+
  geom_point(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+
  geom_line(aes(x=time,  y=Estim,  group=method,  color=method), size=1)+ 
  geom_ribbon(aes(x=time,  ymin=IC_inf,  ymax=IC_sup,  group=method,  fill=method),  alpha=alpha)+ 
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
#' @param color a vector of colors for custom colors, should be the same length
#' as the number of different values of the ic_mat method column.
#' @examples
#' ans <- compute_forall(compute_far.dx_ebm_fit, l_models=c("cnrm", "ipsl"), stat_model=gpd_fit, qthreshold =0.90, p=0.01, R=3)
#' plot_far(ans, axis_trans="al", col=gg_color_hue(2))
#' @export
plot_far <- function(ic_mat, axis_trans, main="", xlim, ylim, col=NULL){
  ticks=c(0, 1/100, 1/10, 1/5, 1/3, 1/2, 1, 2, 3, 5, 10, 100, Inf)
  if(missing(xlim)) xlim <- range(ic_mat$time)
  if(missing(ylim)) ylim <- c(-1.1, 1.1)
  inv <- get(paste(axis_trans, "_inv", sep=""))
  trans <- get(paste(axis_trans, "_trans", sep=""))
  breaks <- trans(ticks) 
  param <- paste(axis_trans,"RR", sep="")
  p1 <- plot_boot_time(ic_mat, param=param, main=main, alpha=0.05) +
  scale_y_continuous(name = "Relative Risk",
                     breaks = breaks,
                     labels = trans_format(inv, function(x) format(x, digits=2)))+
  coord_cartesian(xlim=xlim, ylim=ylim)
  if (!is.null(col)){ 
    p1 <- p1 + scale_color_manual(values=col)
    p1 <- p1 + scale_fill_manual(values=col)
}
  p2 <- p1 + scale_y_continuous(name="FAR",
                                breaks = breaks,
                                labels = trans_format(inv, function(x) format(rrtofar(x), digits=3)))
  ggplot_dual_axis(p1,p2)
}

# plot the relative risk by pannel : each pannel correspond to a diffrent value
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
                     labels = trans_format(inv, function(x) format(x, digits=2)))+
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
