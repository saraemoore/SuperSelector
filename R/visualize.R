################################################################################
# FEATURE SELECTION VISUALIZATION

#' Utility factorizer function
#' 
#' Used by CV SuperLearner feature selection visualization functions
#' 
#' @param feat A vector of values that can be made into a factor
#' @param val A vector by which to order 
#' @param featLabs A named character vector containing factor labels. Can be \code{NULL}.
#' @param rev Should the factor be sorted in ascending order by val? Defaults to \code{TRUE}.
#' @return factorized version of \code{feat}
#' @importFrom dplyr arrange desc summarize group_by
#' @importFrom magrittr `%>%`
featReorderAndFactorize = function(feat, val, featLabs, rev = TRUE) {
    resOrdered = data.frame(feat = feat, val = val) %>%
        group_by(feat) %>%
        dplyr::summarize(val = mean(val))
    resOrdered <- if(rev) {
        arrange(resOrdered, val)
    } else {
        arrange(resOrdered, desc(val))
    }
    resOrdered = as.data.frame(resOrdered)

    featLevels = as.character(feat)
    featLevelsNew = unique(as.character(resOrdered[,"feat"]))

    if(!is.null(featLabs)) {
        featLabs = featLabs[featLevelsNew]
    } else {
        featLabs = featLevelsNew
    }
    factor(featLevels, levels = featLevelsNew, labels = featLabs, ordered = TRUE)
}

#' Short description
#' 
#' Long description
#' 
#' @param res was cvsl.var.imp.all
#' @param feat.labs
#' @param multifold
#' @return A list of ggplots
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by summarize rename filter select
#' @importFrom magrittr `%>%`
#' @importFrom stats median
#' @import ggplot2
cvSLVarImpPlotOld = function(res, feat.labs, multifold = TRUE) {
        # store the ggplots here:
    diag.plots = vector("list", 3)
    # sysfonts::font.add.google("Open Sans", "opensans")

    res = reshape2::melt(res)
    if(multifold) {
        res <- res %>% rename(`cv.fold` = L1, `cv.sl.combo.method` = L2)
    } else {
        # double check this
        res <- res %>% rename(`cv.sl.combo.method` = L1)
    }
    res$cv.sl.combo.method = as.factor(res$cv.sl.combo.method)

    cvsl.var.imp.summ = res %>%
        group_by(cv.sl.combo.method, var) %>%
        dplyr::summarize(mean.val = mean(value))

    feat.labs.df = data.frame(var = names(feat.labs), lab = unname(feat.labs))
    feat.labs.wgtd = merge(cvsl.var.imp.summ %>% filter(cv.sl.combo.method=="wgtd.total") %>% select(var, mean.val),
        feat.labs.df, by="var", all.x=TRUE)
    feat.labs.wgtd = feat.labs.wgtd[order(feat.labs.wgtd$mean.val, decreasing=FALSE),]
    res$lab = factor(as.character(res$var),
        levels=as.character(feat.labs.wgtd$var),
        labels=as.character(feat.labs.wgtd$lab),
        ordered=TRUE)

    # diagnostic plots
    if(multifold){
        diag.plots[[1]] = ggplot(data=subset(res, cv.sl.combo.method=="weighted"), aes(x=value, y=lab, shape = as.factor(cv.fold)))
    }else{
        diag.plots[[1]] = ggplot(data=subset(res, cv.sl.combo.method=="weighted"), aes(x=value, y=lab))
    }
    diag.plots[[1]] = diag.plots[[1]] +
        geom_point() +
        facet_grid(. ~ method, scales = "free_x") +
        # theme_classic(base_family = "opensans") +
        # theme_bw(base_family = "opensans") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 5), legend.position="none") +
        ylab("") +
        ggtitle("Weighted CV SL Variable Importance\n(Algorithm's CV SL coef.) x (prop. of CV folds where feature was selected)")

    if(multifold){
        diag.plots[[3]] = ggplot(data=subset(res, cv.sl.combo.method=="wgtd.total"), aes(x=value, y=lab, shape = as.factor(cv.fold)))
    } else{
        diag.plots[[3]] = ggplot(data=subset(res, cv.sl.combo.method=="wgtd.total"), aes(x=value, y=lab))
    }
    diag.plots[[3]] = diag.plots[[3]] +
        geom_vline(xintercept = c(median(subset(res, cv.sl.combo.method=="wgtd.total")$value),
            mean(subset(res, cv.sl.combo.method=="wgtd.total")$value)), colour="grey80", linetype = "longdash") +
        geom_point() +
        # theme_classic(base_family = "opensans") +
        # theme_bw(base_family = "opensans") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 5), legend.position="none") +
        ylab("") +
        ggtitle("Total weighted CV SL Variable Importance\n(summed over all algorithms)")

    # diag.plots[[2]] = ggplot(data=subset(res, cv.sl.combo.method=="unweighted"), aes(x=value, y=var, group=method, shape=method)) +
    #   geom_point() +
    #   theme_bw() +
    #   theme(axis.text.y = element_text(size = 5)) +
    #   ylab("") +
    #   ggtitle("Unweighted CV SL Variable Importance\n(proportion of CV folds where feature was 'important')") +
    #   scale_shape_manual(values=c(3,4,21,22))
        # scale_shape_manual(values=c(3,4,21,22,24))

    feat.labs.unwgtd = merge(subset(cvsl.var.imp.summ, cv.sl.combo.method=="unwgtd.total")[,c("var","mean.val")], feat.labs.df, by="var", all.x=TRUE)
    feat.labs.unwgtd = feat.labs.unwgtd[order(feat.labs.unwgtd$mean.val, decreasing=FALSE),]
    res$lab = factor(as.character(res$var),
        levels=as.character(feat.labs.unwgtd$var),
        labels=as.character(feat.labs.unwgtd$lab),
        ordered=TRUE)

    if(multifold){
        diag.plots[[2]] = ggplot(data=subset(res, cv.sl.combo.method=="unwgtd.total"), aes(x=value, y=lab, shape = as.factor(cv.fold)))
    }else{
        diag.plots[[2]] = ggplot(data=subset(res, cv.sl.combo.method=="unwgtd.total"), aes(x=value, y=lab))
    }
    diag.plots[[2]] = diag.plots[[2]] +
        geom_vline(xintercept = c(median(subset(res, cv.sl.combo.method=="unwgtd.total")$value),
            mean(subset(res, cv.sl.combo.method=="unwgtd.total")$value)), colour="grey80", linetype = "longdash") +
        geom_point() +
        # theme_bw(base_family = "opensans") +
        # theme_classic(base_family = "opensans") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 5), legend.position="none") +
        ylab("") +
        ggtitle("Total unweighted CV SL Variable Importance\n(summed over all algorithms)")

    if(multifold){
        diag.plots[[1]] = diag.plots[[1]] + scale_shape_manual(values = 48:57)
        diag.plots[[3]] = diag.plots[[3]] + scale_shape_manual(values = 48:57)
        diag.plots[[2]] = diag.plots[[2]] + scale_shape_manual(values = 48:57)
    }

    # showtext::showtext.opts(dpi = 300)
    # showtext::showtext.auto()

    return(diag.plots)
}



    # res = reshape2::melt(res)
    # if(multifold){
    #   colnames(res)[colnames(res)=="L1"] = "cv.fold"
    #   colnames(res)[colnames(res)=="L2"] = "cv.sl.combo.method"
    # }else{
    #   # double check this
    #   colnames(res)[colnames(res)=="L1"] = "cv.sl.combo.method"
    # }
    # res$cv.sl.combo.method = as.factor(res$cv.sl.combo.method)

    # cvsl.var.imp.summ = res %>%
    #   group_by(cv.sl.combo.method, var) %>%
    #   dplyr::summarize(mean.val = mean(value))

    # feat.labs.df = data.frame(var = names(feat.labs), lab = unname(feat.labs))
    # feat.labs.wgtd = merge(subset(cvsl.var.imp.summ, cv.sl.combo.method=="wgtd.total")[,c("var","mean.val")],
    #   feat.labs.df, by="var", all.x=TRUE)
    # feat.labs.wgtd = feat.labs.wgtd[order(feat.labs.wgtd$mean.val, decreasing=FALSE),]
    # res$lab = factor(as.character(res$var),
 #      levels=as.character(feat.labs.wgtd$var),
 #      labels=as.character(feat.labs.wgtd$lab),
 #      ordered=TRUE)
    # if(multifold){
    #   diag.plots[[1]] = ggplot(data=subset(res, cv.sl.combo.method=="weighted"), aes(x=value, y=lab, shape = as.factor(cv.fold)))
    # }else{
    #   diag.plots[[1]] = ggplot(data=subset(res, cv.sl.combo.method=="weighted"), aes(x=value, y=lab))
    # }
    # diag.plots[[1]] = diag.plots[[1]] +
    #   geom_point() +
    #   facet_grid(. ~ method, scales = "free_x") +
    #   # theme_classic(base_family = "opensans") +
    #   # theme_bw(base_family = "opensans") +
    #   theme_bw() +
    #   theme(axis.text.y = element_text(size = 5), legend.position="none") +
    #   ylab("") +
    #   ggtitle("Weighted CV SL Variable Importance\n(Algorithm's CV SL coef.) x (prop. of CV folds where feature was selected)")
    # diag.plots[[1]] = diag.plots[[1]] + scale_shape_manual(values = 48:57)

#' Short description
#' 
#' Long description
#' 
#' @param res
#' @param featCol
#' @param valCol
#' @param valLab
#' @param labelVals
#' @param featLabs
#' @param panelLabs
#' @param varImpTheme
#' @param pointSize
#' @param strokeSize
#' @param maxColor
#' @param title
#' @param subtitle
#' @param panelCol
#' @param shapeCol
#' @param addSummary
#' @return
#' @importFrom scales percent
#' @importFrom stats median
#' @import ggplot2
#' @importFrom ggbeeswarm geom_beeswarm
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
#' res <- cvSuperSelect(Y, X, gaussian(),
#'                      list(SL.library = setNames(list(c("SL.mean", "screen.randomForest.imp"),
#'                                                      c("SL.mean", "screen.earth.backwardprune")),
#'                                                 c("random forest biggest diff mean",
#'                                                   "splines biggest diff mean")),
#'                           method = "method.NNLS",
#'                           selector.library = data.frame(selector = c("cutoff.biggest.diff",
#'                                                                      "cutoff.k"),
#'                                                         k = c(NA, 3),
#'                                                         rowname = c("biggest diff", "top3"),
#'                                                         stringsAsFactors = FALSE) %>%
#'                                              tibble::column_to_rownames(),
#'                           nFolds = 3),
#'                      verbose = TRUE)
#' res_sum <- summarizeScreen(res$summary, groupCols = c("method", "screener"))
#' cvSLVarImpPlot(res_sum)
#' }
cvSLVarImpPlot = function(res,
                          featCol = "term",
                          valCol = "propFoldSel",
                          valLab = "Fraction of CVSL folds in which feature was selected",
                          labelVals = TRUE,
                          featLabs = NULL,
                          panelLabs = NULL,
                          # jitter_vertical = 0.08,
                          varImpTheme = theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                           plot.subtitle = element_text(hjust = 0.5),
                                                           legend.position = "none",
                                                           panel.grid.major.x = element_blank(),
                                                           panel.grid.minor.x = element_blank()),
                          pointSize = 2,
                          strokeSize = 1,
                          maxColor = "grey10",
                          title = "CV SuperLearner Variable Importance",
                          subtitle = NULL,
                          panelCol = NULL,
                          shapeCol = NULL,
                          addSummary = TRUE) {

    if(!valCol %in% colnames(res)) {
        stop("valCol is not a column in res. Check valCol.")
    }

    res[,featCol] = featReorderAndFactorize(res[,featCol], res[,valCol], featLabs)

    if(!is.null(panelCol)&!is.null(panelLabs)) {
        res[,panelCol] <- factor(res[,panelCol], levels = names(panelLabs), labels = panelLabs)
    }
    if(!is.null(shapeCol)) {
        res[,shapeCol] = factor(res[,shapeCol], levels = sort(unique(res[,shapeCol])), ordered = TRUE)
    }

    # library(ggrepel)
    # p <- if(is.null(shapeCol)) {
        # ggplot(data = res, aes_string(x = valCol, y = featCol, label = paste0(".", valCol), shape = shapeCol))
    # } else {
    #   ggplot(data = res, aes_string(x = valCol, y = featCol, label = shapeCol))
    # }

    p <- if(labelVals) {
        # nesting inside aes_string doesn't work...
        res[,paste0(".", valCol)] = as.character(round(res[,valCol], digits = 2))
        ggplot(data = res, aes_string(x = valCol, y = featCol, label = paste0(".", valCol)))
    } else {
        ggplot(data = res, aes_string(x = valCol, y = featCol))
    }

    p = p +
        varImpTheme +
        labs(y = "", title = title, subtitle = subtitle)

    if(!is.null(panelCol)) {
        p = p + facet_grid(paste(".", panelCol, sep = " ~ ")) # , scales = "free_x"
    }

    p <- if(is.null(shapeCol)) {
        p +
            scale_x_continuous(valLab, limits = c(0, 1), labels = scales::percent) +
            # geom_jitter(width = 0, height = jitter_vertical, alpha = 0.6)
            geom_beeswarm(size = pointSize, stroke = strokeSize, color = maxColor, alpha = 0.6, groupOnX = FALSE)
    } else {
        one_code = 49
        p +
            scale_x_continuous(valLab, limits = c(-0.1, 1), labels = scales::percent) +
            # geom_jitter(width = 0, height = jitter_vertical, alpha = 0.6, size = 2) +
            geom_beeswarm(aes_string(shape = shapeCol), size = pointSize, stroke = strokeSize, color = maxColor, alpha = 0.6, groupOnX = FALSE) +
            # geom_point() +
            # geom_text_repel(nudge_x = -0.07, direction = "y", segment.size = 0.2) +
            scale_shape_manual(values = one_code:((one_code - 1) + length(unique(res[,shapeCol]))))
    }
    if(labelVals) {
        p = p + geom_text(vjust = 0, nudge_y = 0.1, check_overlap = TRUE)
    }

    if(addSummary) {
        valSumm = c(median(res[,valCol]), mean(res[,valCol]))
        p = p + geom_vline(xintercept = valSumm, colour = "grey80", linetype = "longdash") +
            geom_text(data = data.frame(x = valSumm, y = 0, label = c("median", "mean")),
                      aes(label = label, x = x, y = y),
                      vjust = 0, hjust = 0, angle = 90)
    }
    return(p)
}

#' Short description
#' 
#' Long description
#' 
#' @param res
#' @param featCol
#' @param featLab
#' @param featLabs
#' @param valCol
#' @param valLab
#' @param catCol
#' @param catLab
#' @param varImpTheme
#' @param pointSize
#' @param strokeSize
#' @param maxColor
#' @param title
#' @param subtitle
#' @param panelCol
#' @param panelLabs
#' @param shapeCol
#' @return
#' @import ggplot2
#' @importFrom ggbeeswarm geom_beeswarm
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
#' res <- cvSuperSelect(Y, X, gaussian(),
#'                      list(SL.library = setNames(list(c("SL.mean", "screen.randomForest.imp"),
#'                                                      c("SL.mean", "screen.earth.backwardprune")),
#'                                                 c("random forest biggest diff mean",
#'                                                   "splines biggest diff mean")),
#'                           method = "method.NNLS",
#'                           selector.library = data.frame(selector = c("cutoff.biggest.diff",
#'                                                                      "cutoff.k"),
#'                                                         k = c(NA, 3),
#'                                                         rowname = c("biggest diff", "top3"),
#'                                                         stringsAsFactors = FALSE) %>%
#'                                              tibble::column_to_rownames(),
#'                           nFolds = 3),
#'                      verbose = TRUE)
#' res_sum <- summarizeScreen(res$summary, groupCols = c("method", "screener"))
#' cvSLVarImpPlot2(res_sum)
#' }
cvSLVarImpPlot2 = function(res,
                           featCol = "term",
                           featLab = "Feature name",
                           featLabs = NULL,
                           valCol = "propFoldSel", #selected
                           valLab = "Proportion of CVSL folds in which feature was selected", #"Feature passed screener"
                           catCol = "screener",
                           catLab = "Screening algorithm", #"Filter and weight combination methods"
                           varImpTheme = theme_classic() + theme(plot.title = element_text(hjust = 0.5),
                                                                 plot.subtitle = element_text(hjust = 0.5),
                                                                 legend.position = "none", #"bottom",
                                                                 axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)), # , size = 8
                           pointSize = 2,
                           strokeSize = 1,
                           maxColor = "grey10",
                           title = "CV SuperLearner Variable Importance",
                           subtitle = NULL,
                           panelCol = NULL,
                           panelLabs = NULL,
                           shapeCol = NULL) {

    res[,featCol] = featReorderAndFactorize(res[,featCol], res[,valCol], featLabs, rev = FALSE)
    res[,catCol] = as.factor(res[,catCol])
    res[,valCol] = as.numeric(res[,valCol])
    if(!is.null(shapeCol)) {
        res[,shapeCol] = factor(res[,shapeCol], levels = sort(unique(res[,shapeCol])), ordered = TRUE)
    }
    if(!is.null(panelCol)) {
        if(is.null(panelLabs)) {
            res[,panelCol] <- as.factor(res[,panelCol])
        } else {
            res[,panelCol] <- factor(res[,panelCol], levels = names(panelLabs), labels = panelLabs)
        }
    }
    # nesting inside aes_string doesn't work...
    # res[,paste0(".", valCol)] = as.character(round(res[,valCol], digits = 2))

    # using hack for alpha gradient legend: https://stackoverflow.com/a/44171042
    # library(scales)
    # ggplot(subset(res3, metafold==1), aes(y = propFoldSel, x = term)) + geom_point() + facet_grid(screener ~ .)

    # p = ggplot(res, aes_string(x = featCol, y = catCol, alpha = valCol, label = paste0(".", valCol))) +
    #   geom_point(aes_string(colour = valCol), alpha = 0) + # create fake gradient legend
    #   geom_point(size = pointSize, colour = maxColor, shape = 18) +
    #   scale_color_gradient(valLab, high = maxColor, low = "white", breaks = c(0, 1), labels = c("Never", "Always")) +
    #   scale_alpha(range = c(0, 1)) +
    #   guides(alpha = FALSE) +
    #   labs(x = featLab, y = catLab, title = title, subtitle = subtitle) +
    #   varImpTheme

    p = ggplot(res, aes_string(x = featCol, y = valCol))
        # geom_point(size = pointSize, stroke = strokeSize, colour = maxColor, shape = 1) +

    # if(length(unique(res[,valCol])) > 2) {
    #   p = p + geom_text(vjust = 0, nudge_y = 0.15, check_overlap = TRUE, size = pointSize*(2/3))
    # }

    p <- if(is.null(shapeCol)) {
        p +
            geom_beeswarm(size = pointSize, stroke = strokeSize, color = maxColor, shape = 1, groupOnX = TRUE) +
            scale_shape(solid = FALSE)
    } else {
        one_code <- 49
        p +
            geom_beeswarm(aes_string(shape = shapeCol), size = pointSize, stroke = strokeSize, groupOnX = TRUE) +
            # scale_colour_brewer(palette = "Paired")
            scale_shape_manual(values = one_code:((one_code - 1) + nlevels(res[,shapeCol])))
    }

    p <- p + labs(x = featLab, y = valLab, title = title, subtitle = subtitle) +
        ylim(0, 1) +
        varImpTheme %+replace% theme(panel.border = element_rect(fill = NA, colour = "black"))

    p <- if(is.null(panelCol)) {
        p + facet_grid(paste(catCol, ".", sep = " ~ "))
    } else {
        p + facet_grid(paste(catCol, panelCol, sep = " ~ "))
    }

    return(p)
}
