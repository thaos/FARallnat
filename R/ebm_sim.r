# used to create function for the EBM simulations
make_ebm_bsamples <- function(boot_ebm){
  # function(mdata, by="year", R, gno2aan=TRUE, ...){
  function(mdata, model,  by="year", R, gno2aan=TRUE){
    # ebm_args <- c(list(...), R=R)
    ebm_args <- list(model=model, R=R)
    bsamples <- do.call(boot_ebm, ebm_args)
    bsamples <- lapply(bsamples, merge, x=mdata, by.y="year", by.x=by)
    if(gno2aan) bsamples <- lapply(bsamples, gno2aan)
    bindexes <- boot_samples(bsamples[[1]], R=R)
    list("bsamples"=bsamples, "bindexes"=bindexes)
  }
}
# function with no scaling factors and the model parameters are taken if available.
# If not, it choses for each bootstrap sample a different sets of parameters
# Example
#' Generate a dataset of climate response to forcings
#' 
#' \code{ebm_bsamples.default} generates the climate responses to external
#' forcings. No random scaling factors are included. The EBM takes the
#' corresponding GCM parameters if available, otherwise it takes an available
#' set of paramters chosent at random
#' @param mdata the data to which the EBM simulations will be added
#' @param model, the name of the GCM for which the climate responses are needed
#' @param by the name of the time variable in mdata. Used to merge mdata with
#' the EBM response.
#' @param R the number of bootstrap simulations
#' @param gno2aan, whether to express the climate response into the ALL, ANT and
#' NAT responses instead of GHG, NAT, and OTHERS(anthropogenic) response
#' @return a list of data.frame where the simulated EBM response are merge to
#' the dataset mdata
#' @examples
#' #Simulate 5 EBM realisations and express the results in ALL, ANT and NAT
#' #responses. In this case, all the EBM simulations are the same since random
#' #scaling factor are not accounted for and because the EBM paramters for the
#' #CNRM are presents
#' data(cnrm)
#' ebm_bsamples.default(mdata=cnrm, R=5, model="cnrm", by="year", gno2aan=TRUE)
#' @export
ebm_bsamples.default <- make_ebm_bsamples(boot_ebm.default)


#' \code{ebm_bsamples.sf_gaussian} generates the climate responses to external
#' forcings with no random scaling factors. it takes the corresponding GCM
#' parameters if available, otherwise it takes the average of the availables
#' parameters
#' @examples
#' # Simulate 5 EBM realisations with no random scaling factors and express the
#' # results in GHG, NAT and OTHERS responses
#' data(cnrm)
#' ebm_bsamples.sf_gaussian(mdata=cnrm, R=5, model="cnrm", by="year", gno2aan=FALSE)
#' @rdname ebm_bsamples.default 
#' @export
ebm_bsamples.sf_gaussian <- make_ebm_bsamples(boot_ebm.cm_default.sf_gaussian)

#' \code{ebm_bsamples.default} generates the climate responses to external
#' forcings. Random gaussian scaling factors on the three categories of forcing
#' (GHG, NAT, and OTHERS ) are included. The EBM takes the
#' corresponding GCM parameters if available, otherwise it takes an available
#' set of paramters chosent at random
#' @examples
#' # Simulate 5 EBM realisations with random scaling factors and express the
#' # results in GHG, NAT and OTHERS responses
#' data(cnrm)
#' ebm_bsamples.sf_gaussian(mdata=cnrm, R=5, model="cnrm", by="year", gno2aan=FALSE)
#' @rdname ebm_bsamples.default 
#' @export
ebm_bsamples.cm_avg.sf_default <- make_ebm_bsamples(boot_ebm.cm_default.sf_gaussian)

# for a dataset, returns the line indexes that will be used for the bootstrap
# samples
boot_samples <- function(data, R=250){
  i_samples <- lapply(1:R, function(x) sample.int(nrow(data), size=nrow(data), replace=TRUE))
  i_samples[[1]] <- 1:nrow(data)
  i_samples
}

# rescale the forcings using random gaussian scaling factors
sf.gaussian <- function(FF_data){
    forcing <- FF_data$FF
    sigma <- FF_data$sigma
    sf <- unlist(lapply(sigma, rnorm, n=1, mean=1))
    for(n in names(sigma))
      forcing[, n] <- forcing[, n] * sf[n]
    forcing
}
# keep the forcings as they are
sf.default <- function(FF_data){FF_data$FF}

# simulates the EBM response for GHG, NAT and OTHERS given a model name and a set
# of forcings
ebm_gno  <- function(model, FF){
    nat <- held_model(FF$nat, model=model)[-1,"T"]
    ghg <- held_model(FF$ghg, model=model)[-1,"T"]
    others <- held_model(FF$others, model=model)[-1,"T"]
    ans <- data.frame("year"=FF$year, "ghg"=ghg, "nat"=nat, "others"=others) 
    # attr(ans, "model") <- model
    ans
}

# creates function to simulate the EBM response given:
# - a function which rescales the forcings (eg. gaussian scaling factors)
# - a function which choses the EBM parameters according to a model name
make_boot_ebm <- function(sf_fun=sf.default, cm_fun=choose_model.default){
  function(model, R){
    data(FF,  envir = environment())
    l_models <- lapply(rep(model, R), cm_fun)
    l_FF <- lapply(seq.int(R), function(r) sf_fun(FF))
    ffdata <- mapply(ebm_gno, model=l_models, FF=l_FF, SIMPLIFY=FALSE)
    ffdata[[1]] <- ebm_gno(l_models[[1]], FF$FF)
    ffdata
  }
}
# no scaling fator, random set of parameters for the EBM if the model is not
# present in set of EBM parameters
boot_ebm.default <- make_boot_ebm(sf_fun=sf.default, cm_fun=choose_model.default)
# random gaussian scaling fator, random set of parameters for the EBM if the model is not
# present in set of EBM parameters
boot_ebm.cm_default.sf_gaussian <- make_boot_ebm(sf_fun=sf.gaussian, cm_fun=choose_model.default)
boot_ebm.cm_avg.sf_default <- make_boot_ebm(sf_fun=sf.default, cm_fun=choose_model.avg)

# express the EBM responses into ALL, NAT and ANT responses
gno2aan <- function(gno_df){
  year <- gno_df$year
  gno_df$ant <- gno_df$ghg + gno_df$others
  gno_df$all <- gno_df$nat + gno_df$ant
  gno_df$ghg <- NULL
  gno_df$others <- NULL
  gno_df
}

# creates function to choose a model if the models has not is own EBM
# parameters
make_choose_model <- function(expr){
  function(model){
    model <- toupper(model)
    l_models <- c("BCC","CCCMA","CNRM","CSIRO","GFDL","INM","IPSL","MIROC","MPIM","MRI","NCC")
    if (model %in% l_models)
      ans <- model
    else
      ans <- eval(expr)
    ans
  }
}
# chooses a set of parameters at radom among the present set of parameters if the
# model does not have its own set of parameters
choose_model.default <- make_choose_model(expression(sample(l_models, 1)))
# chooses a the average set of parameters if the model does not have its own set
# of parameters
choose_model.avg <- make_choose_model(expression("AVG"))


