% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_far_simple.r
\name{make_compute_far_simple}
\alias{make_compute_far_simple}
\title{Creates function to compute the FAR}
\usage{
make_compute_far_simple(ebm_bsamples = ebm_bsamples.default,
  ebm_args = list(model = "random"), decompose_x = dx.gam_allnat,
  dx_args = list())
}
\arguments{
\item{ebm_bsamples}{a function to simulates the EBM responses. it should
return a list of two list : 
\itemize{
\item bsamples : a list of data.frame with the original dataset merged with the
bootstrap EBM simulations
\item bindexes : a list of vectors of indexes indicating the rows that will be
used for the bootstrap
}}

\item{ebm_args}{a list of argument to be used in the ebm_bsamples. It can be
an expression if the variable in the list need to be evaluate within the
compute_far function(Non-Standard-Evaluation NSE)}

\item{decompose_x}{a function to decompose the covariate x into an ALL, an
ANT and a NAT component. It has to takes as argument bsamples and bindexes
which results of the ebm_samples functions}

\item{dx_args}{a list of argument to be used in the ebm_bsamples. It can be}
}
\value{
a function with the following arguments :
\itemize{
\item mdata, a data.frame containg the variables whose names are given by x, y and time
\item y, the name of variable that will be used as the variable of interest y
\item x, the name of variable that will be used as the covariate x 
\item time, the name of variable that will be used as the as the time variable 
\item xp, the threshold used to define the FAR 
\item R, the number of bootstrap samples
\item stat_model the statistical model to explain y in function of x, either
gauss_fit, gev_fit, or gpd_fit from the FARg package
\item ci_p the level of the confidence intervals
\item ... additional arguments if require by the stat_model function
}
}
\description{
Creates function to compute the FAR
}
\examples{
# creates a variante of the computing chain with the following properties 
# EBM simulations: 
# - takes model parameters if available, else takes a set of available
# parameters at random
# - random  scaling factors
# Decompose x:
# -  modèle gam :x_all = beta_nat * nat + s(ant)
# - shift mean(ant) to 0 between 1850 and 1879
compute_far_simple.default <- make_compute_far_simple(ebm_bsamples=ebm_bsamples.default,
                                              ebm_args=expression(list(
                                                                     model=model,
                                                                   mdata=mdata
                                                                 ))
                                      )
library(FARg)
model <- "cnrm"
#load data from the package
data(list=model)
# formating data, e.g passing from temperature to anomalie, keep only hist
# and rcp runs
mdata <- format_data(get(model))
if(model != "obs") mdata <- select_continuous_run(mdata)
ans <- compute_far_simple.default(mdata,
 y="eur_tas", x="gbl_tas", time="year",
  xp=1.6, stat_model=gauss_fit, ci_p=0.9)
}
