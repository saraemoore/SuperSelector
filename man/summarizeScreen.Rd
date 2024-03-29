% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize.R
\name{summarizeScreen}
\alias{summarizeScreen}
\title{Summarize "summary" element of \code{cvSLFeatureSelector}}
\usage{
summarizeScreen(
  res,
  groupCols = c("method", "predictor", "screener"),
  discrete = FALSE
)
}
\arguments{
\item{res}{The \code{data.frame} "summary" element of a
\link{cvSLFeatureSelector} result}

\item{groupCols}{Column names of \code{res} (as a character vector) to group
by before summarizing.}

\item{discrete}{If \code{TRUE}, instead of summarizing, simply return rows in \code{res} where
\code{discrete} is \code{TRUE}. Defaults to \code{FALSE}.}
}
\value{
A \code{data.frame}
}
\description{
Summarize "summary" element of \code{cvSLFeatureSelector}
}
\examples{
\dontrun{
# based on example in SuperLearner package
dat <- sim_sl_data(n_obs = 100, rnd_seed = 1)
res <- cvSLFeatureSelector(dat \%>\% pull(Y), dat \%>\% select(-c(ID, Y)), family = gaussian(),
                           method = "method.NNLS",
                           SL.library = setNames(list(c("SL.mean", "screen.randomForest.imp"),
                                                      c("SL.mean", "screen.earth.backwardprune")),
                                                 c("random forest biggest diff mean",
                                                   "splines biggest diff mean")),
                           selector.library = data.frame(selector = c("cutoff.biggest.diff",
                                                                     "cutoff.k"),
                                                         k = c(NA, 3),
                                                         rowname = c("biggest diff", "top3"),
                                                         stringsAsFactors = FALSE) \%>\%
                                              tibble::column_to_rownames(),
                           nFolds = 3,
                           verbose = TRUE)
summarizeScreen(res$summary, groupCols = c("method", "screener"))
summarizeScreen(res$summary, groupCols = "method")
}
}
