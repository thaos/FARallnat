% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constrained_far.r
\name{constrained_far}
\alias{constrained_far}
\title{Creates function to co compute the FAR}
\usage{
constrained_far(far_o, far_m, model, xp = 1.6, ci_p = 0.9)
}
\arguments{
\item{far_o}{result from a function compute_far}

\item{far_m}{result from a function compute_far. It should have the same
stastical model for y as far_o. the fitted parameters with the exception of
the intercepts will constrained the statistical model that will be refitted
for far_o.}
}
\value{
the far_o object where the the fitted model far_o$mfit has been
constrain by the results of far_m
}
\description{
\code{constrained_far} compute the far for wich the non-stationnary
components from model fit_o are constrained by the statistical model fit_m.
i.e. fit_o must have the same fitted parameters has fit_m with the exception
of the intercept parameters.
}
\examples{
far_cnrm <- compute_far.default("cnrm", y="eur_tas", x="gbl_tas", time="year")
far_obs <- compute_far.default("obs", y="eur_tas", x="gbl_tas", time="year")
far_cstr <- constrained_far(far_obs, far_cnrm, model="cnrm")
}
