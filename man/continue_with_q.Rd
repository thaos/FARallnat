% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.r
\name{continue_with_q}
\alias{continue_with_q}
\title{ans <- compute_forall(compute_far.default, stat_model=gauss_fit, p=0.01, R=3)
ans_cstr <- continue_with_constrained_far(ans) 
summary_plot(ans_cstr)
Compute the quantile corresponding to a probability of exceedance p}
\usage{
continue_with_q(ans, model, p)
}
\arguments{
\item{ans}{results of a compute_far function}

\item{model}{the name of the model for which compute_far has been run}

\item{p}{the probability used to define the quantile}
}
\value{
a list with the same objects as a compute_far function with in
addition an object allq which contains all the bootstrap estimates of q
}
\description{
\code{continue_with_q} computes the quantile corresponding to a probability
of exceedance p from the results of a compute_far function
}
\examples{
# compute the FAR for the CNRM model using a  gam decomposition
# and a gaussian fit with only three bootstrap samples
ans <- compute_far.default(model="cnrm", y="eur_tas", x="gbl_tas", time="year", xp=1.6, stat_model=gauss_fit, ci_p=0.9)
ans <- continue_with_q(ans, model="cnrm", p=0.01)
}
