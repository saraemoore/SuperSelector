#' Short description
#' 
#' Long description
#' 
#' @param res The \code{summary} element of the result of CV SuperLearner feature selection.
#' @param groupCols
#' @param collapseCols
#' @param discrete If \code{TRUE}, instead of summarizing, simply return rows in \code{res} where
#' \code{discrete} is \code{TRUE}. Defaults to \code{FALSE}. 
#' @return A \code{data.frame}
#' @importFrom dplyr filter group_by group_by_at summarize arrange desc
#' @importFrom magrittr `%>%`
#' @export 
#' @examples
#' \dontrun{
#' # based on example in SuperLearner package
#' set.seed(1)
#' n <- 100
#' p <- 20
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' X <- data.frame(X)
#' Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + rnorm(n)
#' res <- cvSLFeatureSelector(Y, X, family = gaussian(),
#'                            method = "method.NNLS",
#'                            SL.library = setNames(list(c("SL.mean", "screen.randomForest.imp"),
#'                                                       c("SL.mean", "screen.earth.backwardprune")),
#'                                                  c("random forest biggest diff mean",
#'                                                    "splines biggest diff mean")),
#'                            selector.library = data.frame(selector = c("cutoff.biggest.diff",
#'                                                                      "cutoff.k"),
#'                                                          k = c(NA, 3),
#'                                                          rowname = c("biggest diff", "top3"),
#'                                                          stringsAsFactors = FALSE) %>%
#'                                               tibble::column_to_rownames(),
#'                            nFolds = 3,
#'                            verbose = TRUE)
#' summarizeScreen(res$summary, groupCols = c("method", "screener"))
#' summarizeScreen(res$summary, groupCols = c("method", "screener"), collapseCols = "screener")
#' }
summarizeScreen = function(res, groupCols = c("method", "predictor", "screener"), collapseCols = NULL, discrete = FALSE) {

    if(discrete) {
        # basically do nothing
        return(res %>% filter(discrete))
    }

    # summarize allResByFold$selected by method+screener+term (proportion of folds term/variable was selected)
    for(g in groupCols) {
        if((g %in% colnames(res))&!(g %in% collapseCols)) {
            res = res %>% group_by_at(g, .add = TRUE)
        }
    }
    res = res %>%
            group_by(term, .add = TRUE)

    res <- if("predictor" %in% colnames(res)) {
        res %>%
            dplyr::summarize(propFoldSel = mean(selected),
                      # numFoldSel = sum(selected),
                      wtdPropFoldSel = mean(estimate*selected)) %>%
            arrange(desc(wtdPropFoldSel))
                      # wtdNumFoldSel = sum(estimate*selected)))
    } else {
         res %>%
            dplyr::summarize(propFoldSel = mean(selected)) %>%
            arrange(desc(propFoldSel))
                      # numFoldSel = sum(selected))
    }

    res %>%
        as.data.frame

    # TODO: make term an ordered factor according to numFoldSel?
    # TODO: make method a factor?
}
