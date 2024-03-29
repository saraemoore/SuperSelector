% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{make_dag_toy_data}
\alias{make_dag_toy_data}
\title{Create DAG for a toy example}
\usage{
make_dag_toy_data()
}
\value{
A \code{DAG} object as defined by the \code{simcausal} package.
}
\description{
Create DAG for a toy example of 3 continuous covariates and 1 binary outcome where the outcome
(\code{Y}) is unrelated to the first covariate (normally-distributed \code{W1}) but is related
to the second and third covariates (lognormally-distributed \code{W2} and \code{W3}, which are
also related to one another). Note that this function requires the \code{simcausal} package,
which may need to be installed from github like so:
\code{remotes::install_github('osofr/simcausal', build_vignettes = FALSE)}
}
\examples{
\dontrun{
d <- make_dag_toy_data()
}
}
