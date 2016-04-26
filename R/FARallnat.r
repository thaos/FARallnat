#' FARallnat  : A Toolbox package for the Analysis of the FAR on the European
#' Summer temperature.
#'
#' The aims is explain the distribution changes of a variable y with respect to
#' variations of a variable x
#'
#' The FAR analysis is decompose in three big steps:
#'
#' @section Simulations with an EBM of the climate responses: 
#'
#' The simulations are done with respect to three categories of
#' forcings : anthropogenic(ANT), natural(NAT), and others(OTHERS)
#' Three variantes are implemented :
#' \itemize{
#' \item \code{emb_samples.default}
#' \item \code{emb_samples.sf_gaussian}
#' \item \code{emb_samples.cm_avg.sf_default}
#'}
#'
#' @section Decomposition of the variable x into ALL,  ANT and NAT:
#'
#' Most of the times the EBM simulations help for the decomposition.
#' several variantes are implemented :
#' \itemize{
#' \item \code{dx.gam_allnat}
#' \item \code{dx.gam_gno}
#' \item \code{dx.lm_gno}
#' \item \code{dx.raw}
#' \item \code{dx.ebm_fit}
#'} 
#'
#' @section Fit a statistical model to y:
#' 
#' The model explains variations in the distributions of y with with respect to
#' the covariate x. From this fit quantities such as the FAR and
#' the probilities of exceedance can be computed. For the fit, we use three
#' functions of the FARg package
#' \itemize{
#' \item \code{gauss_fit}
#' \item \code{gev_fit}
#' \item \code{gpd_fit}
#'} 
#' The foo functions ...
#'
#' @section Combining the three steps:
#'
#' Functions that combine those three steps can be created with the
#' function \code{make_compute_far}. Once a compute_far function is created, it
#' can be applied to all the available datasets of CMIP-5 simulations using the
#' function \code{compute_forall} and a first visualization can be done with
#' \code{summary_plot}
#'
#' @docType package
#' @name FARallnat 
NULL
