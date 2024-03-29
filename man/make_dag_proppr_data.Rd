% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{make_dag_proppr_data}
\alias{make_dag_proppr_data}
\title{Create DAG for a PROPPR example}
\usage{
make_dag_proppr_data(outcome_name = "Y")
}
\arguments{
\item{outcome_name}{Character string containing the name of the outcome
variable. Defaults to \code{"Y"} if not specified.}
}
\value{
A \code{DAG} object as defined by the \code{simcausal} package.
}
\description{
Create DAG for a PROPPR example of 9 covariates and 1 binary outcome where the outcome
is unrelated to the first four covariates (normally-distributed \code{W1:4}) but is related
to the remaining five covariates. Note that this function requires the \code{simcausal} package,
which may need to be installed from github like so:
\code{remotes::install_github('osofr/simcausal', build_vignettes = FALSE)}
}
\examples{
\dontrun{
d <- make_dag_proppr_data()
}
}
