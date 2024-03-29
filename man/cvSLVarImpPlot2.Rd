% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualize.R
\name{cvSLVarImpPlot2}
\alias{cvSLVarImpPlot2}
\title{An alternative visualization of the summarized result from \code{cvSLFeatureSelector}}
\usage{
cvSLVarImpPlot2(
  res,
  featCol = "term",
  featLab = "Feature name",
  featLabs = NULL,
  valCol = "propFoldSel",
  valLab = "Proportion of CVSL folds in which feature was selected",
  catCol = "screener",
  catLab = "Screening algorithm",
  varImpTheme = theme_classic(base_size = 10) + theme(plot.title = element_text(hjust =
    0.5), plot.subtitle = element_text(hjust = 0.5), legend.position = "none",
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)),
  pointSize = 2,
  strokeSize = 1,
  maxColor = "grey10",
  title = "CV SuperLearner Variable Importance",
  subtitle = NULL,
  panelCol = NULL,
  panelLabs = NULL,
  shapeCol = NULL,
  y_breaks = waiver(),
  tikz_output = FALSE
)
}
\arguments{
\item{res}{A \code{data.frame} containing at least the columns specified in
\code{featCol}, \code{valCol}, and \code{catCol} and, if specified,
\code{panelCol} and \code{shapeCol}.}

\item{featCol}{A character string corresponding to a column name of
\code{res} containing the feature/predictor variable name. Defaults to
"term".}

\item{featLab}{An optional character string containing the horizontal axis
label. Defaults to "Feature name".}

\item{featLabs}{An optional named character vector containing feature labels.}

\item{valCol}{A character string corresponding to a column name of
\code{res} containing the value to be plotted on the vertical axis.
Defaults to "propFoldSel".}

\item{valLab}{A character string containing the vertical axis label.
Defaults to "Proportion of CVSL folds in which feature was selected".}

\item{catCol}{A character string corresponding to a column name of
\code{res} containing the variable to use for vertical plot paneling.
Defaults to "screener".}

\item{catLab}{A character string containing the vertical plot paneling
label. Defaults to "Screening algorithm".}

\item{varImpTheme}{A ggplot theme for the plot.}

\item{pointSize}{The point size for the plot. Defaults to \code{2}.}

\item{strokeSize}{The line thickness for the plot. Defaults to \code{1}.}

\item{maxColor}{The darkest color of the scatterplot points. Defaults to
\code{"grey10"}.}

\item{title}{An optional title for the plot. Defaults to
"CV SuperLearner Variable Importance".}

\item{subtitle}{An optional subtitle for the plot. If not provided, no
subtitle will appear on the plot.}

\item{panelCol}{An optional character string corresponding to a column name
of \code{res} containing the variable to use for horizontal plot paneling. If
omitted, the resulting plot will only be paneled vertically.}

\item{panelLabs}{An optional named character vector containing labels for
the vertical panels.}

\item{shapeCol}{An optional character string corresponding to a column name
of \code{res} containing the variable to use for scatterplot point shapes.
If omitted, a single default shape will be used. If set, typically set to
"fold".}

\item{y_breaks}{An optional vector of breaks to use on the vertical axis.}

\item{tikz_output}{A boolean indicating whether percent symbols should be
escaped. Defaults to \code{FALSE}.}
}
\value{
A \code{ggplot}
}
\description{
An alternative visualization of the summarized result from \code{cvSLFeatureSelector}
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
                           nFolds = 5,
                           verbose = TRUE)
res_sum <- summarizeScreen(res$summary, groupCols = c("method", "screener"))
cvSLVarImpPlot2(res_sum)
}
}
