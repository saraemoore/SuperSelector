
SuperSelector
=============

<!-- badges: start -->

[![R-CMD-check](https://github.com/saraemoore/SuperSelector/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/saraemoore/SuperSelector/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Overview
--------

The `SuperSelector` package extends Cross-Validated
[SuperLearner](https://github.com/ecpolley/SuperLearner)
(`SuperLearner::CV.SuperLearner()`) for ensembled feature selection.
This package also contains utility functions for simulating input data
for, and visualizing outputs from, its feature selection algorithm.

Installation
------------

The SuperSelector package is currently only available via GitHub. To
install:

``` r
remotes::install_github("saraemoore/SuperSelector")
```

Examples
--------

### Simple

From `?sim_toy_data`:

> Simulate data for a 3 continuous covariate, 1 binary outcome example
> where the outcome (Y) is unrelated to the first covariate
> (normally-distributed W1) but is related to the second and third
> covariates (lognormally-distributed W2 and W3, which are also related
> to one another).

``` r
library(SuperSelector)
library(SuperLearner)  # for SL.mean, method.NNloglik
library(FSelector)     # for cutoff.biggest.diff

dat <- sim_toy_data(n_obs = 1000, rnd_seed = 54321)
res <- cvSLFeatureSelector(dat$Y,
                           dat[,-which(colnames(dat) %in% c("ID", "Y"))],
                           family = binomial())
```

``` r
knitr::kable(res$whichVariable[,"keep_names", drop = FALSE])
```

|                          | keep\_names |
|:-------------------------|:------------|
| cutoff.biggest.diff\_kNA | W2          |

``` r
res_sum <- summarizeScreen(res$summary,
                           groupCols = "method")
p1 <- cvSLVarImpPlot(res_sum)
p1
```

<img src="man/figures/README-simple_example_plot1-1.png" width="100%" />

``` r
res_sum2 <- summarizeScreen(res$summary,
                            groupCols = c("method", "screener"))
p2 <- cvSLVarImpPlot2(res_sum2)
p2
```

<img src="man/figures/README-simple_example_plot2-1.png" width="100%" />

### Advanced

From `?sim_proppr_data`:

> Simulate data for a 9 covariate, 1 binary outcome example where the
> outcome (Y) is unrelated to the first four covariates
> (normally-distributed W1:4) but is related to the remaining five
> covariates (some of which are also related to one another).

``` r
library(SuperSelector)
library(SuperLearner)  # for SL.mean, SL.glm, method.NNloglik
library(FSelector)     # for cutoff.biggest.diff, cutoff.k
# remotes::install_github("saraemoore/SLScreenExtra")
library(SLScreenExtra) # for screen.wgtd.lasso, screen.randomForest.imp, screen.earth.backwardprune

dat <- sim_proppr_data(n_obs = 500, rnd_seed = 54321)
dat_x <- dplyr::bind_cols(dat[, -which(colnames(dat) %in% c("GCS", "ID", "Y"))],
                          factor_to_indicator("GCS", dat))

libraryCVSLFeatSel <- list(
    `lasso mean` = c("SL.mean", "screen.wgtd.lasso"),
    `random forest biggest diff mean` = c("SL.mean", "screen.randomForest.imp"),
    `splines biggest diff mean` = c("SL.mean", "screen.earth.backwardprune"),
    `lasso glm` = c("SL.glm", "screen.wgtd.lasso"),
    `random forest biggest diff glm` = c("SL.glm", "screen.randomForest.imp"),
    `splines biggest diff glm` = c("SL.glm", "screen.earth.backwardprune")
)

libraryMetaFeatSel <- data.frame(selector = c("cutoff.biggest.diff", "cutoff.k", "cutoff.k"),
                                k = c(NA, 3, 6),
                                stringsAsFactors = FALSE)
rownames(libraryMetaFeatSel) <- c("biggest diff", "top3", "top6")

res <- cvSLFeatureSelector(as.numeric(dat$Y),
                           dat_x,
                           family = binomial(),
                           SL.library = libraryCVSLFeatSel,
                           selector.library = libraryMetaFeatSel,
                           nFolds = 5)
```

``` r
knitr::kable(res$whichVariable[,"keep_names", drop = FALSE])
```

|              | keep\_names                                         |
|:-------------|:----------------------------------------------------|
| biggest diff | GCS.15                                              |
| top3         | GCS.15 , PrehospOCY , GCS.4\_to\_14                 |
| top6         | GCS.15 , PrehospOCY , GCS.4\_to\_14, Age , W1 , BMI |

``` r
res_sum <- summarizeScreen(res$summary,
                           groupCols = "method")
p1 <- cvSLVarImpPlot(res_sum)
p1
```

<img src="man/figures/README-advanced_example_plot1-1.png" width="100%" />

``` r
res_sum2 <- summarizeScreen(res$summary,
                            groupCols = c("method", "screener"))
p2 <- cvSLVarImpPlot2(res_sum2)
p2
```

<img src="man/figures/README-advanced_example_plot2-1.png" width="100%" />
