% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{sim_proppr_data}
\alias{sim_proppr_data}
\title{Simulate data for a PROPPR example}
\usage{
sim_proppr_data(
  n_obs = 500,
  rnd_seed = NULL,
  outcome_name = "Y",
  verbose = FALSE
)
}
\arguments{
\item{n_obs}{Number of observations to simulate.}

\item{rnd_seed}{Random seed to be used when simulating.}

\item{outcome_name}{Name of outcome column (only impacts DAG output)}

\item{verbose}{Print diagnostic messages? Defaults to FALSE}
}
\value{
A \code{data.frame} with columns \code{ID}, \code{W1:4}, \code{GCS}, \code{Age},
\code{PrehospOCY}, \code{BMI}, \code{SysBP}, and \code{Y}.
}
\description{
Simulate data for a 9 covariate, 1 binary outcome example where the outcome
(\code{Y}) is unrelated to the first four covariates (normally-distributed \code{W1:4}) but is related
to the remaining five covariates (some of which are
also related to one another). Note that this function requires the \code{simcausal} package,
which may need to be installed from github like so:
\code{remotes::install_github('osofr/simcausal', build_vignettes = FALSE)}
}
\examples{
\dontrun{
dat <- sim_proppr_data(n_obs = 200, rnd_seed = 620)
testobs <- sim_proppr_data(n_obs = 1, rnd_seed = 801)
}
}
