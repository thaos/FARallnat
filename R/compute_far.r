# -------------------------------------------------------------------------------------
#' Creates function to co compute the FAR 
#'
#' \code{make_compute_far} creates a function to compute the FAR by choosing an
#' EBM simulation function and a method to decompose the covariate x into an ANT and a NAT
# component
#'
#' @param ebm_bsamples a function to simulates the EBM responses. it should
#' return a list of two list : 
#' \itemize{
#' \item bsamples : a list of data.frame with the original dataset merged with the
#' bootstrap EBM simulations
#' \item bindexes : a list of vectors of indexes indicating the rows that will be
#' used for the bootstrap
#' }
#' @param ebm_args a list of argument to be used in the ebm_bsamples. It can be
#' an expression if the variable in the list need to be evaluate within the
#' compute_far function(Non-Standard-Evaluation NSE)
#' @param decompose_x a function to decompose the covariate x into an ALL, an
#' ANT and a NAT component. It has to takes as argument bsamples and bindexes
#' which results of the ebm_samples functions
#' @param dx_args a list of argument to be used in the ebm_bsamples. It can be
#' @return a function with the following arguments :
#' \itemize{
#' \item model, the name of the model to load the data from
#' \item y, the name of variable that will be used as the variable of interest y
#' \item x, the name of variable that will be used as the covariate x 
#' \item time, the name of variable that will be used as the as the time variable 
#' \item xp, the threshold used to define the FAR 
#' \item R, the number of bootstrap samples
#' \item stat_model the statistical model to explain y in function of x, either
#' gauss_fit, gev_fit, or gpd_fit from the FARg package
#' \item ci_p the level of the confidence intervals
#' \item ... additional arguments if require by the stat_model function
#' }
#' @examples
#' # creates a variante of the computing chain with the following properties 
#' # EBM simulations: 
#' # - takes model parameters if available, else takes a set of available
#' # parameters at random
#' # - random  scaling factors
#' # Decompose x:
#' # -  modèle gam :x_all = beta_nat * nat + s(ant)
#' # - shift mean(ant) to 0 between 1850 and 1879
#' compute_far.default <- make_compute_far(ebm_bsamples=ebm_bsamples.default,
#'                                               ebm_args=expression(list(
#'                                                                      model=model,
#'                                                                    mdata=mdata
#'                                                                  ))
#'                                       )
#' ans <- compute_far.default(model="cnrm", y="eur_tas", x="gbl_tas", time="year", xp=1.6, stat_model=gauss_fit, ci_p=0.9)
#' @export
make_compute_far <- function(ebm_bsamples=ebm_bsamples.default,
                                 ebm_args=list(model="cnrm"),
                                 decompose_x=dx.gam_allnat,
                                 dx_args=list()){ 
  # computing chain for the FAR
  function(model, y="y", x="x", time="time", xp=1.6 , R=3, stat_model=gauss_fit, ci_p=0.9, ...){
    # load data from the package
    data(list=model)
    # formating data, e.g passing from temperature to anomalie, keep only hist
    # and rcp runs
    mdata <- format_data(get(model))
    # keep only run that posses rcp simulations
    if(model != "obs") mdata <- select_continuous_run(mdata)
    # -------------------------------------------------------------------------------
    # add EBM simulated response to the datasets
    ebm_bsamples <- match.fun(ebm_bsamples)
    ebm_args <- eval(ebm_args)
    ebm_args <- c(ebm_args, R=R)
    bsamples <- do.call(ebm_bsamples, ebm_args)
    # --------------------------------------------------------------------------------
    # decompose the covariate x into x_all = x_ant + x_nat 
    decompose_x <- match.fun(decompose_x)
    dx_args <- eval(dx_args)
    dx_args <- c(dx_args, y=y, x=x, time=time, bsamples=list(bsamples$bsamples), list(bindexes=bsamples$bindexes))
    dx <- do.call(decompose_x, dx_args)
    #---------------------------------------------------------------------------------
    # fit the statical model : stat_model=gauss_fit, gpd_fit or gev_fit
    # and compute probability of exceending the threshold xp
    mfit <- fit_and_boot_allnat(dx$l_x_an, stat_model, ...)
    allp <- boot_allp_onperiod(mfit, dx$l_x_an_origin, l_time=sort(unique(dx$l_x_an_origin[[1]]$time)), xp=xp)
    # computing the Relative Risk(RR) for each bootstrap sample
    far <- add_param(allp, operation=p_all/p_nat, name="RR")
    # if RR=NA remplace by either +inf or 0
    far["RR",,] %<>% imput_aurel(.) 
    # add rescaled far (notably, atan(log(RR)) the FAR plot 
    far %<>%  add_param(., operation=al_trans(RR), name="alRR") %>% 
    add_param(., operation=el_trans(RR), name="elRR")
    # compute  confidence intervals for each parameters
    ic_far <- get_ic_onperiod(far, method_name=model, ci_p=ci_p)  
    ic_far$param <- as.character(ic_far$param)
    # ----------------------------------------------------------------------------------
    # return a list with notably all bootsrap decomposition dx and all
    # statistical fits for y.
    list("y"=y, "x"=x, "time"=time, "data"=mdata, "dx"=dx, "mfit"= mfit, "allp"=allp, "ic_far"=ic_far)
  } 
}

# -------------------------------------------------------------------------------------
#' List of already created compute_far functions
#' \itemize{
#' \item compute_far.default
#' EBM simulations: takes model parameters if available, otherwise takes an
#' available set of parameters at randoms; no scaling factors 
#'Decompose x: gam model :x_all = beta_nat * nat + s(ant); shift mean(ant) to 0 between 1850 and 1879
#'}
#' @rdname make_compute_far 
#' @export
compute_far.default <- make_compute_far(ebm_bsamples=ebm_bsamples.default,
                                                ebm_args=expression(list(
                                                                         model=model,
                                                                         mdata=mdata
                                                                         ))
                                                )
#' \itemize{
#' \item compute_far.sf_gaussian
#' EBM simulations: takes model parameters if available, otherwise takes an
#' available set of parameters at randoms; random scaling factors 
#'Decompose x: gam model :x_all = beta_nat * nat + s(ant); shift mean(ant) to 0 between 1850 and 1879
#'}
#' @rdname make_compute_far 
#' @export
compute_far.sf_gaussian <- make_compute_far(ebm_bsamples=ebm_bsamples.sf_gaussian,
                                                    ebm_args=expression(list(
                                                                             model=model,
                                                                             mdata=mdata
                                                                             ))
                                                    )

#' \itemize{
#' \item compute_far.sf_gaussian.dx_raw
#' EBM simulations: takes model parameters if available, otherwise takes an
#' available set of parameters at randoms; random scaling factors 
#'Decompose x: keeps EBM simulations as they are 
#'}
#' @rdname make_compute_far 
#' @export
compute_far.sf_gaussian.dx_raw <- make_compute_far(ebm_bsamples=ebm_bsamples.sf_gaussian,
                                                           ebm_args=expression(list(
                                                                                    model=model,
                                                                                    mdata=mdata
                                                                                    )),
                                                           decompose_x=dx.raw
                                                           )


#' \itemize{
#' \item compute_far.sf_gaussian.dx_lm_gno
#' EBM simulations: takes model parameters if available, otherwise takes an
#' available set of parameters at randoms; random scaling factors ; responses ghg, nat, and others instead of ant and nat
#'Decompose x: linear model: x_all = beta_nat * nat + beta_ghg * ghg
#'+ beta_others * others; shift mean(ghg) and mean(others) to 0 between 1850 and 1879
#'}
#' @rdname make_compute_far 
#' @export
compute_far.sf_gaussian.dx_lm_gno <- make_compute_far(ebm_bsamples=ebm_bsamples.sf_gaussian,
                                                           ebm_args=expression(list(
                                                                                    model=model,
                                                                                    mdata=mdata,
                                                                                    gno2aan=FALSE
                                                                                    )),
                                                           decompose_x=dx.lm_gno
                                                           )

#' \itemize{
#' \item compute_far.sf_gaussian.dx_gam_gno
#' EBM simulations: takes model parameters if available, otherwise takes an
#' available set of parameters at randoms; random scaling factors ; responses ghg, nat, and others instead of ant and nat
#'Decompose x: gam model : x_all = beta_nat * nat + beta_ghg * ghg
#'+ s(others); shift mean(ghg) and mean(others) to 0 between 1850 and 1879
#'}
#' @rdname make_compute_far 
#' @export
compute_far.sf_gaussian.dx_gam_gno <- make_compute_far(ebm_bsamples=ebm_bsamples.sf_gaussian,
                                                           ebm_args=expression(list(
                                                                                    model=model,
                                                                                    mdata=mdata,
                                                                                    gno2aan=FALSE
                                                                                    )),
                                                           decompose_x=dx.gam_gno
                                                           )

#' \itemize{
#' \item compute_far.dx_ebm_fit
#' EBM simulations: takes model parameters if available, otherwise takes an
#' available set of parameters at randoms; no scaling factors
#' Decompose x: EBM_fit (ordinary least square)
#'}
#' @rdname make_compute_far 
#' @export
compute_far.dx_ebm_fit <- make_compute_far(ebm_bsamples=ebm_bsamples.default,
                                                           ebm_args=expression(list(
                                                                                    model=model,
                                                                                    mdata=mdata,
                                                                                    gno2aan=FALSE
                                                                                    )),
                                                           decompose_x=dx.ebm_fit
                                                           )
