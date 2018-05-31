# -------------------------------------------------------------------------------------
#' function to compute the FAR
#'
#' \code{compute_far_simple} Function to compute the FAR from annual time-series
# -------------------------------------------------------------------------------------
#' @param mdata, a data.frame containing the variables whose names are given by x, y and time
#' @param y, the name of variable that will be used as the variable of interest y
#' @param x, the name of variable that will be used as the covariate x 
#' @param time, the name of variable that will be used as the as the time variable 
#' @param xp, the threshold used to define the FAR 
#' @param R, the number of bootstrap samples
#' @param stat_model the statistical model to explain y in function of x, either
#' gauss_fit, gev_fit, or gpd_fit from the FARg package
#' @param ci_p the level of the confidence intervals
#' @param ... additional arguments if required by the stat_model function
#' @examples
#' library(FARg)
#' #load data from the package
#' # formating data, e.g passing from temperature to anomalie, keep only hist
#' # and rcp runs
#' mdata <- format_data(cnrm)
#' # keep only run that posses rcp simulations
#' mdata <- select_continuous_run(mdata)
#' ans <- compute_far_simple(mdata,
#'  y="eur_tas", x="gbl_tas", time="year",
#'   xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' @export
compute_far_simple <- function(mdata, y="y", x="x", time="time", xp=1.6 , R=3, stat_model=gauss_fit, ci_p=0.9, ...){
    # -------------------------------------------------------------------------------
    # add EBM simulated response to the datasets
    bsamples <- ebm_bsamples.default(mdata = mdata, model = "random", R = R)
    # --------------------------------------------------------------------------------
    # decompose the covariate x into x_all = x_ant + x_nat 
    dx <- dx.gam_allnat(bsamples$bsamples, bsamples$bindexes, y = y, x = x, time = time)
    #---------------------------------------------------------------------------------
    # fit the statical model : stat_model=gauss_fit, gpd_fit or gev_fit
    # and compute probability of exceending the threshold xp
    mfit <- fit_and_boot_allnat(dx$l_x_an, stat_model, ...)
    allp <- boot_allp_onperiod(mfit, dx$l_x_an_origin, l_time=sort(unique(dx$l_x_an_origin[[1]]$time)), xp=xp)
    # computing the Relative Risk(RR) for each bootstrap sample
    far <- add_param(allp, operation=p_all/p_nat, name="RR")
    # if RR=NA remplace by either +inf or 0
    far["RR",,] <- imput_aurel(far["RR",,]) 
    # add rescaled far (notably, atan(log(RR)) the FAR plot 
    far  <- add_param(far, operation=al_trans(RR), name="alRR") 
    far <-  add_param(far, operation=el_trans(RR), name="elRR")
    # compute  confidence intervals for each parameters
    ic_far <- get_ic_onperiod(far, ci_p=ci_p)  
    ic_far$param <- as.character(ic_far$param)
    # ----------------------------------------------------------------------------------
    # return a list with notably all bootsrap decomposition dx and all
    # statistical fits for y.
    list("y"=y, "x"=x, "time"=time, "data"=mdata, "dx"=dx, "mfit"= mfit, "allp"=allp, "ic_far"=ic_far)
}
