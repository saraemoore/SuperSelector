################################################################################
# SCREENING/FEATURE SELECTION via CV.SuperLearner

#' A selector which uses an absolute cutoff
#' 
#' In the style of \code{\link[FSelector]{cutoff.k}}.
#' 
#' @param attrs A \code{data.frame} containing ranks in the first column and attribute names as
#' row names. See \code{\link[FSelector]{cutoff.k}}.
#' @param k A numeric threshold between the min and max of the first column of \code{attrs}.
#' @return Rownames of the first column of \code{attrs} that exceed \code{k}.
cutoff.val = function(attrs, k) {
# TODO: should the values first be normalized before comparing?
#       is this even valid? should the cutoff instead be some percentile (median?) or mean?
    if (dim(attrs)[1] == 0)
        return(character(0))
    if (k < min(attrs[,1]))
        warning("k too small. No attrs will be removed.")
    if (k >= max(attrs[,1])) {
        stop("k too large.")
    }

    attrs = attrs[order(attrs[,1], decreasing = TRUE),,drop = FALSE]
    rownames(attrs)[attrs[,1] > k]
}

#' Short description
#' 
#' Long description
#' 
#' @param attrs
#' @param selector
#' @param k
#' @param env Environment in which \code{selector} can be found. Defaults to \code{parent.frame()}.
#' @return The result of applying \code{selector} with argument \code{k} to \code{attrs}.
featureSelector = function(attrs, selector, k, env = parent.frame()) {
    selector_f <- get(selector, envir = env)
    # c(list()) trick so k disappears if NULL
    do.call(selector_f, c(list(attrs = attrs), k = k))
}

#' Coerce a value to numeric or die trying
#' 
#' Coerce a value to numeric or die trying
#' 
#' @param x A value that should be coerced to numeric
#' @param errorFunction A function which elicits a useful error message if \code{x} cannot be
#' coerced to numeric
#' @return \code{x} as numeric, if it was coercible; otherwise, an error is produced.
#' @export
#' @examples
#' \dontrun{
#' errfxn = function(txt = "") {
#'   stop(paste("x must be NA, NULL, or coercible to numeric.", txt))
#' }
#' x <- "5"
#' x <- makeLibArgNumeric(x, errfxn)
#' }
makeLibArgNumeric = function(x, errorFunction) {
    if(is.na(x))
        x = NULL
    if(is.character(x)) {
        x = tryCatch(as.numeric(x),
                     error = errorFunction,
                     warning = errorFunction)
    }
    if(!(is.numeric(x)|is.null(x)))
        errorFunction()

    return(x)
}

#' Short description
#' 
#' Long description
#' 
#' @param selector
#' @param k
#' @param x
#' @param xNames
#' @param env environment; passthru to featureSelector. defaults to parent.frame()
#' @param verbose logical; if TRUE, print retained features. defaults to FALSE
#' @importFrom stats setNames
#' @return
metaFeatSel = function(selector, k, x, xNames, env = parent.frame(), verbose = FALSE) {

    attrs = as.data.frame(x)
    rownames(attrs) = xNames

    if(!is.character(selector)) {
        "selector must be a character string"
    }

    kStop = function(txt = "") {
        stop(paste("k must be NA, NULL, or coercible to numeric.", txt))
    }
    k = makeLibArgNumeric(k, kStop)

    keep <- featureSelector(attrs, selector, k, env)

    if(verbose) {
        cat("FEATURES KEPT:\n", paste(keep, collapse=", "), "\n")
    }

    return(setNames(xNames %in% keep, xNames))
}

############################################################
# wrap the below in a method called extractScreen.CV.SuperLearner
# also create an S3 generic extractScreen
# and a default method extractScreen.default
# this method could be implemented for all screening algorithms (including those in SLScreenExtra)

#' Short description
#' 
#' Long description
#' 
#' @param df
#' @param sel
#' @param propCol
#' @param verbose logical; passthru to metaFeatSel
#' @importFrom purrr map2 map
#' @importFrom dplyr mutate
#' @importFrom magrittr `%>%`
#' @return
selectFeaturesBySelector <- function(df, sel, propCol, verbose) {
    # as.character because it may come in as a factor if stringsAsFactors isn't FALSE
    sel %>%
        mutate(keep = map2(as.character(selector), k, metaFeatSel, df[, propCol], df[, "term"], verbose = verbose)) %>%
        mutate(keep_names = map(keep, function(x) names(which(x))))
}

#' Short description
#' 
#' Long description
#' 
#' @param df
#' @param sel
#' @param propCol
#' @param verbose logical; passthru to selectFeaturesBySelector
#' @importFrom purrr map
#' @importFrom dplyr bind_rows
#' @importFrom magrittr `%>%`
#' @return
selectFeaturesByMethodAndSelector <- function(df, sel, propCol, verbose) {
    df %>%
        split(.$method) %>%
        map(~ selectFeaturesBySelector(.x, sel, propCol, verbose = verbose)) %>%
        bind_rows(.id = "method")
}

#' Short description
#' 
#' Long description
#' 
#' @param x
#' @param selector.library
#' @param weighted
#' @param verbose logical; passthru to selectFeaturesByMethodAndSelector
#' @return
#' @import broom
#' @import SLScreenExtra
#' @importFrom dplyr bind_rows
extractScreen.CV.SuperLearner = function(x, selector.library, weighted, verbose) {

    # alg <- ifelse(weighted, "both", "screening")
    allResByMethod = lapply(c("both", "screening"), function(alg)
        sapply(x, tidy, algorithm = alg, simplify = FALSE))
    names(allResByMethod) <- c("weighted", "unweighted")

    allRes <- lapply(allResByMethod, dplyr::bind_rows, .id = "method")

    # collapseAcross <- if(weighted) {
    #   # if weighted, additionally collapse across prediction algorithm.
    #   c("predictor", "screener")
    # } else {
    #   # summarize proportion of times each term was chosen
    #   # collapse across fold and screening algorithm.
    #   "screener"
    # }
    summByFeature <- mapply(summarizeScreen,
                            allRes, collapseCols = list(c("predictor", "screener"), "screener"),
                            SIMPLIFY = FALSE)
    # for output purposes, better to have this as a list of data.frames than a single data.frame.
    # summByFeature = split(summByFeature %>% select(-method), summByFeature$method)

    # propCol <- ifelse(weighted, "wtdPropFoldSel", "propFoldSel")

    if(is.null(rownames(selector.library))|rownames(selector.library)==seq(nrow(selector.library))) {
        rownames(selector.library) = apply(selector.library, 1, paste, collapse = "_k")
    }

    cvSLkeepCols <- mapply(function(df, propCol)
                           selectFeaturesByMethodAndSelector(df, selector.library, propCol, verbose = verbose),
                           summByFeature, c("wtdPropFoldSel", "propFoldSel"),
                           SIMPLIFY = FALSE)
    # names(cvSLkeepCols) <- c("weighted", "unweighted")

    return(list(whichVariable = cvSLkeepCols[[c("unweighted", "weighted")[weighted+1]]],
                summary = allRes[[c("unweighted", "weighted")[weighted+1]]]))
}

############################################################

#' Short description
#' 
#' Long description
#' 
#' @param x a vector
#' @return logical
is_not_homogenous <- function(x) {
    return(length(unique(x)) > 1)
}

#' Short description
#' 
#' Long description
#' 
#' @param wv_result The \code{data.frame} stored in the \code{whichVariable} element of the list
#' returned from \code{extractScreen.CV.SuperLearner}. This \code{data.frame} is expected to
#' contain columns named \code{keep} and \code{method}. 
#' @param feature_names Ordered names of predictor variable(s).
#' @return A \code{data.frame} with indicators reordered by original \code{features} column order
#' @importFrom dplyr mutate rename
#' @importFrom magrittr `%>%`
#' @importFrom purrr map map_chr
sortWhichVariable <- function(wv_result, feature_names) {
    wv_result %>%
        mutate(keep = map(keep, function(i) i[feature_names])) %>%
        mutate(keep_bin = map_chr(keep, function(x) paste(as.numeric(x), collapse = ""))) %>%
        rename(combo_method = method) # prevent conflicts later
}

#' Feature selection via CV.SuperLearner
#' 
#' Feature selection via CV.SuperLearner
#' 
#' @param Y Outcome (numeric vector). See \code{\link[SuperLearner]{CV.SuperLearner}}.
#' @param X Predictor variable(s) (data.frame or matrix). See
#' \code{\link[SuperLearner]{CV.SuperLearner}}.
#' @param family Error distribution to be used in the model:
#' \code{\link[stats]{gaussian}} or \code{\link[stats]{binomial}}.
#' See \code{\link[SuperLearner]{CV.SuperLearner}}.
#' @param obsWeights Optional numeric vector of observation weights. See
#' \code{\link[SuperLearner]{CV.SuperLearner}}.
#' @param id Cluster identification variable. See
#' \code{\link[SuperLearner]{CV.SuperLearner}}.
#' @param method A list of method(s) by which \code{\link[SuperLearner]{CV.SuperLearner}} should
#' estimate coefficients to combine algorithms.
#' @param SL.library A list of character vectors of length 2, each containing a screener algorithm
#' and a prediction algorithm. See \code{\link[SuperLearner]{SuperLearner}}.
#' @param selector.library As SL.library, either a vector of length 2 or a list of such vectors. The first element of each vector should be a string naming a selector function (typically from the \code{\link[FSelector]{FSelector}} package, with a prefix of \code{cutoff.}). The second element of the vector should be the required second argument to the named function (usually named \code{k}), if one exists. If there is no second argument, this second element can be set to \code{NULL} or, equivalently, omitted altogether.
#' @param nFolds numeric of length 1 or 2. If length(nFolds)==1, the value provided will be used as the number of cross-validation folds for both the outer (\code{\link[SuperLearner]{CV.SuperLearner}}) and inner (SuperLearner) cross-validations. If length(nFolds)==2, the first element or element with name "outer" will be used as the number of folds for the outer (\code{\link[SuperLearner]{CV.SuperLearner}}) cross-validation, and the second element or element with name "inner" will be used as the number of folds for the inner (SuperLearner) cross-validation.
#' @param trimLogit Only applicable when using the \code{NNloglik} method. See
#' \code{\link[SuperLearner]{SuperLearner.control}}.
#' @param stratifyCV Only applicable when \code{Y} is binary. If \code{TRUE},
#' \code{\link[SuperLearner]{CV.SuperLearner}} will stratify CV splits by \code{Y}. See
#' \code{\link[SuperLearner]{SuperLearner.CV.control}}.
#' @param shuffle If \code{TRUE}, rows of \code{X} will be shuffled before being split. See
#' \code{\link[SuperLearner]{SuperLearner.CV.control}}.
#' @param validRows A list containing pre-specified rows for the CV splits. See
#' \code{\link[SuperLearner]{SuperLearner.CV.control}}.
#' @param weighted Should the weights estimated by the method(s) be used to weight the feature
#' selection? Passed to \link{extractScreen.CV.SuperLearner}.
#' @param verbose Print diagnostic messages? Defaults to FALSE
#' @param label An optional named character vector of length 1. If specified, the value will be added as a column (where the column name is set to \code{names(label)}) in the \code{data.frame} stored in the \code{summary} element of the \code{cvslFull} element of the returned list. One example of when this might be useful is when this function is called from within a cross-validation fold. Then, \code{label} might be set to, for example, \code{c(metafold = fold$v)}.
#' @param ... Passed through to \code{\link[SuperLearner]{CV.SuperLearner}}
#' @return A named list containing the results of the \code{\link[SuperLearner]{CV.SuperLearner}}
#' feature selection. Will contain elements \code{whichVariable} (a \code{data.frame}),
#' \code{summary} (a \code{data.frame}), and cvslFull (a \code{list} containing one result of class
#' \code{CV.SuperLearner} for each `method` supplied).
#' @importFrom SuperLearner CV.SuperLearner recombineCVSL
#' @importFrom dplyr select
#' @importFrom tidyselect all_of
#' @importFrom magrittr `%>%`
#' @importFrom stats binomial gaussian
#' @importFrom tibble add_column
#' @export
#' @examples
#' \dontrun{
#' # remotes::install_github('osofr/simcausal', build_vignettes = FALSE)
#' dat <- sim_toy_data(n_obs = 200, rnd_seed = 620)
#' res <- cvSLFeatureSelector(dat %>% pull(Y), dat %>% select(-c(ID, Y)), family = binomial(),
#'                            method = "method.NNloglik",
#'                            SL.library = setNames(list(c("SL.mean", "screen.randomForest.imp"),
#'                                                       c("SL.mean", "screen.earth.backwardprune")),
#'                                                  c("random forest biggest diff mean",
#'                                                    "splines biggest diff mean")),
#'                            selector.library = data.frame(selector = c("cutoff.biggest.diff",
#'                                                                       "cutoff.k"),
#'                                                          k = c(NA, 2),
#'                                                          rowname = c("biggest diff", "top2"),
#'                                                          stringsAsFactors = FALSE) %>%
#'                                               tibble::column_to_rownames(),
#'                            nFolds = 3,
#'                            verbose = TRUE)
#' 
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
#' }
cvSLFeatureSelector = function(Y, X, family = binomial(), obsWeights = NULL, id = NULL, method = "method.NNloglik",
    SL.library = list(c("SL.mean", "screen.corP"), c("SL.mean", "screen.glmnet"), c("SL.mean", "screen.randomForest")),
    selector.library = data.frame(selector = "cutoff.biggest.diff", k = NA, stringsAsFactors = FALSE),
    nFolds = c(outer = 10, inner = 10),
    trimLogit = 0.001, stratifyCV = (family$family=="binomial"), shuffle = TRUE, validRows = NULL, weighted = FALSE,
    verbose = FALSE, label = NULL, ...) {

    if(length(nFolds) < 1) {
        nFolds = c(outer = 10, inner = 10)
    } else if(length(nFolds) == 1) {
        nFolds = rep(nFolds, 2)
    } else if(length(nFolds) > 2) {
        stop("nFolds cannot contain more than 2 elements.")
    }
    if(is.null(names(nFolds))) {
        names(nFolds) = c("outer", "inner")
    }

    # need to be able to supply more than one value for weighted.

    if(!weighted) {
        if(length(method) > 1) {
            warning("Multiple library algorithm combination methods supplied to argument 'method'\n",
                    " but 'weighted' is FALSE. Unweighted results will not change across methods.\n",
                    " Therefore, only the first element of 'method' will be used.")
        }
        # method doesn't matter unless looking at coefficient estimates
        method = method[[1]]
    }

    # remove columns from X if they're homogeneous
    not_homogenous_cols <- X %>% select(where(is_not_homogenous)) %>% colnames()
    if(length(not_homogenous_cols) < ncol(X)) {
        warning("The following homogeneous features were dropped in cvSLFeatureSelector():\n",
                paste(setdiff(colnames(X), not_homogenous_cols), collapse = ", "))
        X <- X %>% select(all_of(not_homogenous_cols))
    }

    cvSLres = CV.SuperLearner(
        Y = Y, X = X, family = family,
        SL.library = unname(SL.library),
        method = method[[1]], id = id, verbose = verbose,
        control = list(saveFitLibrary = TRUE, trimLogit = trimLogit), # do we really need saveFitLibrary here?
        cvControl = list(V = nFolds["outer"], stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows),
        innerCvControl = list(list(V = nFolds["inner"], stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)),
        obsWeights = obsWeights, saveAll = TRUE, ...)

    if(weighted&length(method) > 1) {
        cvSLres = c(list(cvSLres),
                    lapply(method[-1], function(i)
                           recombineCVSL(cvSLres,
                           method = i, verbose = verbose, saveAll = TRUE,
                           parallel = "seq")))
    } else {
        cvSLres = list(cvSLres)
    }
    names(cvSLres) = sub("method.", "", method)

    # don't need this for any code below as of now... keep in case needed in future
    # if(is.null(names(SL.library))) {
        # names(SL.library) = cvSLres[[1]]$libraryNames
    # }

    screenRes = extractScreen.CV.SuperLearner(cvSLres,
                                              selector.library = selector.library,
                                              weighted = weighted,
                                              verbose = verbose)

    # reorder indicators by original X column order
    screenRes$whichVariable = sortWhichVariable(screenRes$whichVariable, colnames(X))

    if(!is.null(label)&!is.null(names(label))) {
        # optionally: label with, for example, fold number
        screenRes$summary <- screenRes$summary %>% add_column(!!!label)
    }

    return(c(screenRes, cvslFull = list(cvSLres)))
}

################################################################################

# need to retain method, selector, k -- as three list-columns
# only keeping each unique combo of X columns
# BUT need to keep up with method, selector, k for each so we can expand back out after prediction step

#' Short description
#' 
#' Long description
#' 
#' @param res returned from cvSLFeatureSelector()
#' @return
#' @importFrom dplyr group_by summarize_at left_join select
#' @importFrom magrittr `%>%`
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
#' groupBySelectionSet(res)
#' }
groupBySelectionSet <- function(res) {
    wv <- res$whichVariable
    wv %>%
        group_by(keep_bin) %>%
        summarize_at(c("combo_method", "selector", "k"), list) %>%
        left_join(wv %>% select(keep, keep_names, keep_bin) %>% unique(), by = "keep_bin")
}

################################################################################

#' Wraps \code{cvSLFeatureSelector} for use as a \code{SuperLearner} screening algorithm
#' 
#' Long description
#' 
#' @param Y
#' @param X
#' @param family
#' @param obsWeights
#' @param id
#' @param ...
#' @return
#' @export
#' @examples
#' \dontrun{
#' screen.CV.SuperLearner(Y, X, family = binomial(), obsWeights = NULL, id = NULL, verbose = TRUE)
#' }
# TODO: should this instead go in extra_sl_libalgs.R?
screen.CV.SuperLearner = function(Y, X, family, obsWeights, id, ...) {
    # SL.library, selector.library, etc. can be passed through ...
    res = cvSLFeatureSelector(Y, X, family, obsWeights, id, ...)
    # only keep first set of variables selected just in case
    return(res$whichVariable[[1, "keep"]])
}