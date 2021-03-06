# load the data for one model, turn it into 1961-1990 anomalies
load_data <- function(model){
  tas <- data(list=model,  envir = environment())
  tas <- get(model)
  tas$year <- as.numeric(tas$year)
  tas_subset <- subset(tas, year >= 1961 & year <= 1990 & type != "nat")
  eur_bias <- tapply(as.numeric(tas_subset$eur_tas), tas_subset$run, mean)
  gbl_bias <- tapply(as.numeric(tas_subset$gbl_tas), tas_subset$run, mean)
  uruns <- unique(tas_subset$run)
  for(i in seq_along(uruns)){
    r <- uruns[i]
    tas[tas$run==r, "eur_tas"] <- as.numeric(tas[ tas$run==r, "eur_tas"]) - eur_bias[i] 
    tas[tas$run==r, "gbl_tas"] <- as.numeric(tas[ tas$run==r, "gbl_tas"]) - gbl_bias[i] 
  }
  tas$hnat <- tas$type=="nat"
  # tas$ant[tas$hnat] <- 0
  tas
}

#' Function to prepare the tas datasets contained in the package
#'
#' \code{format_data} turns tas data into 1961-1991 anomalies
#' @param tas : one of tas dataset included in the package
#' @export
format_data <- function(tas){
  tas$year <- as.numeric(tas$year)
  tas_subset <- subset(tas, year >= 1961 & year <= 1990 & type != "nat")
  eur_bias <- tapply(as.numeric(tas_subset$eur_tas), tas_subset$run, mean)
  gbl_bias <- tapply(as.numeric(tas_subset$gbl_tas), tas_subset$run, mean)
  uruns <- unique(tas_subset$run)
  for(i in seq_along(uruns)){
    r <- uruns[i]
    tas[tas$run==r, "eur_tas"] <- as.numeric(tas[ tas$run==r, "eur_tas"]) - eur_bias[i] 
    tas[tas$run==r, "gbl_tas"] <- as.numeric(tas[ tas$run==r, "gbl_tas"]) - gbl_bias[i] 
  }
  tas$hnat <- tas$type=="nat"
  # tas$ant[tas$hnat] <- 0
  tas
}

#' Function to prepare the tas datasets contained in the package
#'
#' \code{select_continuous_run} keeps only climate runs with historical and rcp
#' @param tas : one of tas dataset included in the package
#' @export
select_continuous_run <- function(tas){
  tas <- tas[tas$hnat==0, ]
  runs <- unique(tas[, c("run", "type")])
  runs <- tapply(tas$type, tas$run, function(x)length(unique(x)) >= 2)
  runs <- as.numeric(names(runs)[runs])
  tas[tas$run %in% runs, ]
}

