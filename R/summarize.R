#' Summarize "summary" element of \code{cvSLFeatureSelector}
#' 
#' @param res The \code{data.frame} "summary" element of a
#' \link{cvSLFeatureSelector} result
#' @param groupCols Column names of \code{res} (as a character vector) to group
#' by before summarizing.
#' @param discrete If \code{TRUE}, instead of summarizing, simply return rows in \code{res} where
#' \code{discrete} is \code{TRUE}. Defaults to \code{FALSE}. 
#' @return A \code{data.frame}
#' @importFrom dplyr filter group_by group_by_at summarize arrange desc
#' @importFrom magrittr `%>%`
#' @importFrom rlang .data
#' @export 
#' @examples
#' \dontrun{
#' # based on example in SuperLearner package
#' dat <- sim_sl_data(n_obs = 100, rnd_seed = 1)
#' res <- cvSLFeatureSelector(dat %>% pull(Y), dat %>% select(-c(ID, Y)), family = gaussian(),
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
#' summarizeScreen(res$summary, groupCols = "method")
#' }
summarizeScreen = function(res,
                           groupCols = c("method", "predictor", "screener"),
                           discrete = FALSE) {

    if(discrete) {
        # basically do nothing
        return(res %>% filter(discrete))
    }

    # summarize allResByFold$selected by method+screener+term (proportion of folds term/variable was selected)
    for(g in groupCols) {
        if((g %in% colnames(res))) {
            res = res %>% group_by_at(g, .add = TRUE)
        }
    }
    res = res %>%
            group_by(.data$term, .add = TRUE)

    res <- if("predictor" %in% colnames(res)) {
        res %>%
            dplyr::summarize(propFoldSel = mean(.data$selected),
                      # numFoldSel = sum(.data$selected),
                      wtdPropFoldSel = mean(.data$estimate*.data$selected)) %>%
            arrange(desc(.data$wtdPropFoldSel))
                      # wtdNumFoldSel = sum(.data$estimate*.data$selected)))
    } else {
         res %>%
            dplyr::summarize(propFoldSel = mean(.data$selected)) %>%
            arrange(desc(.data$propFoldSel))
                      # numFoldSel = sum(.data$selected))
    }

    res %>%
        as.data.frame

    # TODO: make term an ordered factor according to numFoldSel?
    # TODO: make method a factor?
}
