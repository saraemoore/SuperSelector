#' A custom smoother for a \code{GGally::ggpairs()} plot
#' 
#' @param data data set using
#' @param mapping aesthetics being used
#' @param ... other arguments to add to \code{\link[ggplot2]{geom_point}}
#' @param method parameters supplied to \code{\link[ggplot2]{geom_smooth}}
#' @return A \code{ggplot}
#' @importFrom ggplot2 ggplot geom_point geom_smooth
ggally_custom_smooth <- function (data, mapping, ..., method = "lm") {
	# initial covariate space:
	# override GGally version with se = FALSE

    p <- ggplot(data = data, mapping)
    p <- p + geom_point(...)
    if (!is.null(mapping$color) || !is.null(mapping$colour)) {
        p <- p + geom_smooth(method = method, se = FALSE)
    }
    else {
        p <- p + geom_smooth(method = method, colour = I("black"), se = FALSE)
    }
    p$type <- "continuous"
    p$subType <- "smooth"
    return(p)
}

#' A custom density for a \code{GGally::ggpairs()} plot
#' 
#' @param data data set using
#' @param mapping aesthetics being used
#' @param ... other arguments to add to \code{\link[ggplot2]{geom_density}}
#' @return A \code{ggplot}
#' @importFrom ggplot2 ggplot scale_y_continuous geom_density
ggally_custom_density <- function(data, mapping, ...){
    p <- ggplot(data, mapping) + scale_y_continuous()
    p <- p + geom_density(...) # fill = NA
    return(p)
}

#' Plot SPLOM for data with binary outcome
#' 
#' Plot SPLOM for \code{data.frame} \code{dat}. 
#' 
#' @param dat A \code{data.frame} with outcome column \code{Y} and ID column \code{ID}.
#' @param outcome_levels Unique allowed values in the outcome column (\code{Y}).
#' @param outcome_colors Vector of colors with a length equal to \code{length(outcome_levels)}.
#' @return SPLOM visualization (ggplot)
#' @importFrom GGally ggpairs wrap grab_legend
#' @importFrom ggplot2 ggplot aes_string geom_hex scale_fill_manual scale_colour_manual
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr mutate select
#' @export
#' @examples
#' \dontrun{
#' dat <- sim_toy_data(n_obs = 200, rnd_seed = 620)
#' p <- plot_splom_sim_data(dat)
#' }
plot_splom_sim_data <- function(dat, outcome_levels = 0:1, outcome_colors = RColorBrewer::brewer.pal(3, "RdBu")[c(1,3)]) {
	if(!all(unique(dat$Y)%in%outcome_levels)) {
		stop("'outcome_levels' does not match dat$Y. Check value passed as 'outcome_levels.'")
	}
	if(length(unique(dat$Y)) > 8) {
		stop("plot_splom_sim_data() not implemented for continuous Y.")
	}

	dat <- dat %>%
		select(-ID) %>%
		mutate(Y = factor(Y, levels = outcome_levels, labels = as.character(outcome_levels)))

	x_names <- dat %>% select(-Y) %>% colnames()
	prettyLegend <- ggplot(dat, aes_string(x = x_names[1], y = x_names[2], fill = "Y")) +
	    geom_hex() +
	    scale_fill_manual("Y", values = outcome_colors)

	p <- ggpairs(dat,
				 columns = which(colnames(dat)!="Y"),
				 mapping = aes(color = Y),
				 title = "Simulation summary",
				 upper = "blank",
				 # upper = list(continuous = wrap("density", alpha = 0.5),
				 # 			 combo = "box",
				 # 			 discrete = wrap("facetbar", alpha = 0.5)),
				 diag = list(continuous = ggally_custom_density),
				 lower = list(continuous = wrap(ggally_custom_smooth, alpha = 0.5, shape = 21, stroke = 1.5),
	                      	  combo = wrap("dot", alpha = 0.5)),
				 legend = grab_legend(prettyLegend))

	for (row in seq_len(p$nrow)) {
	    for (col in seq_len(p$ncol)) {
	      if(!is.null(p[row, col])) {
	        p[row, col] <- p[row, col] +
	            scale_colour_manual("Y", values = outcome_colors)
	      }
	    }
	}

	return(p)
}

#' A heatmap displaying observation-level features and outcome
#' 
#' @param dat A \code{data.frame} with outcome column \code{Y} and ID column \code{ID}.
#' @return A \code{ggplot}
#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_gradient2 labs scale_x_discrete scale_y_discrete theme_bw theme element_blank element_text
#' @importFrom scales brewer_pal
#' @importFrom dplyr select pull
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @export
#' @examples
#' \dontrun{
#' dat <- sim_toy_data(n_obs = 200, rnd_seed = 620)
#' p <- plot_heatmap_sim_data(dat)
#' }
plot_heatmap_sim_data <- function(dat) {
	Y <- dat %>% pull("Y")
	if(length(unique(Y)) > 8) {
		stop("plot_heatmap_sim_data() not implemented for continuous Y.")
	}

	X <- dat %>% select(!all_of(c("ID", "Y")))
	outcome_colors <- scales::brewer_pal("qual")(length(unique(Y)))
	
	# heatmap of expression: green = suppressed, red = overexpressed, black = average
	ggplot() +
	 	geom_tile(data = reshape2::melt(as.matrix(X)), 
	 			  aes(x = as.factor(.data[["Var1"]]), y = .data[["Var2"]], fill = .data[["value"]]),
	 			  color = "white") +
		scale_fill_gradient2("Value", low = "#01665e", mid = "#f5f5f5", high = "#b35806") +
		labs(x = "Observation", y = "Variable") +
	    scale_x_discrete(expand = c(0, 0)) + 
	    scale_y_discrete(expand = c(0, 0)) +
	    theme_bw(base_size = 9) +
	    theme(panel.grid.major.x = element_blank(),
	          panel.grid.major.y = element_blank(),
	          axis.ticks = element_blank(),
	          axis.text.x = element_text(size = 4, face = "bold", angle = 90, hjust = 0.3,
								         # 'Vectorized input to `element_text()` is not officially supported.''
								         color = outcome_colors[Y+1]))
}

################################################################################

#' Simulate data as in a SuperLearner example
#' 
#' Simulate data for a 20 continuous (normally-distributed) covariate, 1 continuous outcome example
#' where the outcome (\code{Y}) is only related to the first three covariates (\code{X1}, \code{X2},
#' and \code{X3}) and no covariate is related to any other covariate.
#' 
#' @param n_obs Number of observations to simulate.
#' @param rnd_seed Random seed to be used when simulating.
#' @return A \code{data.frame} with columns \code{ID}, \code{X1:20}, and \code{Y}.
#' @importFrom magrittr `%>%`
#' @importFrom dplyr bind_cols
#' @importFrom stats rnorm
#' @export
#' @examples
#' \dontrun{
#' dat <- sim_sl_data(n_obs = 200, rnd_seed = 620)
#' testobs <- sim_sl_data(n_obs = 1, rnd_seed = 801)
#' }
sim_sl_data <- function(n_obs = 100, rnd_seed = NULL) {
	set.seed(rnd_seed)

	# based on example in SuperLearner package
	p_feat <- 20
	X <- rnorm(n_obs*p_feat) %>%
		matrix(nrow = n_obs, ncol = p_feat) %>%
		data.frame()
	epsilon <- rnorm(n_obs)
	Y <- X[, 1] + sqrt(abs(X[, 2] * X[, 3])) + X[, 2] - X[, 3] + epsilon
	bind_cols(data.frame(ID = seq_along(Y)),
			  X,
			  data.frame(Y = Y))
}

################################################################################

#' Create DAG for a PROPPR example
#' 
#' Create DAG for a PROPPR example of 9 covariates and 1 binary outcome where the outcome
#' is unrelated to the first four covariates (normally-distributed \code{W1:4}) but is related
#' to the remaining five covariates. Note that this function requires the \code{simcausal} package,
#' which may need to be installed from github like so:
#' \code{remotes::install_github('osofr/simcausal', build_vignettes = FALSE)} 
#' 
#' @param outcome_name Character string containing the name of the outcome
#' variable. Defaults to \code{"Y"} if not specified.
#' @return A \code{DAG} object as defined by the \code{simcausal} package.
#' @importFrom simcausal DAG.empty node set.DAG
#' @importFrom stats plogis
#' @examples
#' \dontrun{
#' d <- make_dag_proppr_data()
#' }
make_dag_proppr_data <- function(outcome_name = "Y") {
	# estimated from empirical distribution of PROPPR data
	d <- DAG.empty() +
	  	node("W1", distr="rnorm", mean = -2, sd = 1) +
	  	node("W2", distr="rnorm", mean = -1, sd = 1) +
	  	node("W3", distr="rnorm", mean = 1, sd = 1) +
	  	node("W4", distr="rnorm", mean = 2, sd = 1) +
		# simplify to 3 categories: gcs=3, gcs>3&gcs<15, and gcs=15.
	  	node("GCS", distr="rcat.factor", probs = c(0.32, 0.28, 0.40)) +
	  	node("Age", distr="rlnorm", meanlog = 3.56, sdlog = 0.44) +
	  	node("PrehospOCY", distr="rexp", rate = 0.54) +
	  	node("BMI", distr="rlnorm", meanlog = 0.08 * log(Age) + 3.01, sdlog = 0.21) +
	  	node("SysBP", distr="rnbinom", mu = exp(-0.002 * Age + 0.016 * PrehospOCY + 4.7), size = 10.9) +
		node(outcome_name, distr="rbern", prob = plogis(
			-2.68 + ifelse(GCS == 1, 1.34, ifelse(GCS == 3, -0.45, 0)) + 0.00012 * SysBP - 0.013 * BMI + 0.015 * Age - 0.16 * PrehospOCY))
	set.DAG(d)
}

#' Plot DAG for PROPPR example
#' 
#' @return DAG visualization
#' @importFrom simcausal plotDAG
#' @export
#' @examples
#' \dontrun{
#' p <- plot_dag_proppr_data()
#' }
plot_dag_proppr_data <- function() {
	d <- make_dag_proppr_data()
	plotDAG(d, xjitter = 0.3, yjitter = 0.01,
			edge_attrs = list(width = 0.5, arrow.width = 0.7, arrow.size = 0.5),
			vertex_attrs = list(size = 12, label.cex = 0.8, shape="none", 
			label.color = "grey20", label.family = "sans"))
}

#' Simulate data for a PROPPR example
#' 
#' Simulate data for a 9 covariate, 1 binary outcome example where the outcome
#' (\code{Y}) is unrelated to the first four covariates (normally-distributed \code{W1:4}) but is related
#' to the remaining five covariates (some of which are
#' also related to one another). Note that this function requires the \code{simcausal} package,
#' which may need to be installed from github like so:
#' \code{remotes::install_github('osofr/simcausal', build_vignettes = FALSE)} 
#' 
#' @param n_obs Number of observations to simulate.
#' @param rnd_seed Random seed to be used when simulating.
#' @param outcome_name Name of outcome column (only impacts DAG output)
#' @param verbose Print diagnostic messages? Defaults to FALSE
#' @return A \code{data.frame} with columns \code{ID}, \code{W1:4}, \code{GCS}, \code{Age},
#' \code{PrehospOCY}, \code{BMI}, \code{SysBP}, and \code{Y}.
#' @importFrom simcausal simobs
#' @export
#' @examples
#' \dontrun{
#' dat <- sim_proppr_data(n_obs = 200, rnd_seed = 620)
#' testobs <- sim_proppr_data(n_obs = 1, rnd_seed = 801)
#' }
sim_proppr_data <- function(n_obs = 500, rnd_seed = NULL, outcome_name = "Y", verbose = FALSE) {
	d <- make_dag_proppr_data(outcome_name = outcome_name)
	dat <- simobs(d, n = n_obs, rndseed = rnd_seed, verbose = verbose)
	dat$GCS <- factor(as.numeric(dat$GCS), levels=1:3, labels=c("3", "4 to 14", "15"))
	return(dat)
}

################################################################################

#' Create DAG for a toy example
#' 
#' Create DAG for a toy example of 3 continuous covariates and 1 binary outcome where the outcome
#' (\code{Y}) is unrelated to the first covariate (normally-distributed \code{W1}) but is related
#' to the second and third covariates (lognormally-distributed \code{W2} and \code{W3}, which are
#' also related to one another). Note that this function requires the \code{simcausal} package,
#' which may need to be installed from github like so:
#' \code{remotes::install_github('osofr/simcausal', build_vignettes = FALSE)} 
#' 
#' @return A \code{DAG} object as defined by the \code{simcausal} package.
#' @importFrom simcausal DAG.empty node set.DAG
#' @importFrom stats plogis
#' @examples
#' \dontrun{
#' d <- make_dag_toy_data()
#' }
make_dag_toy_data <- function() {
	d <- DAG.empty() +
	    node("W1", distr = "rnorm", mean = 1, sd = 1) + # was random noise
	    node("W2", distr = "rlnorm", meanlog = 3.56, sdlog = 0.44) + # was age
	    node("W3", distr = "rlnorm", meanlog = 0.08 * log(W2) + 3.01, sdlog = 0.21) + # was BMI
	    node("Y", distr = "rbern", prob = plogis(-2.35 - 0.013 * W3 + 0.015 * W2))
	set.DAG(d)
}

#' Plot DAG for a toy example
#'
#' @return DAG visualization
#' @importFrom simcausal plotDAG
#' @export
#' @examples
#' \dontrun{
#' p <- plot_dag_toy_data()
#' }
plot_dag_toy_data <- function() {
	d <- make_dag_toy_data()
	plotDAG(d, xjitter = 0.5, yjitter = 0.8,
    		edge_attrs = list(width = 0.3, arrow.width = 0.3, arrow.size = 0.5),
    		vertex_attrs = list(size = 12, label.cex = 1.25, shape="none",
        	label.color = "grey20"))
}

#' Simulate data for a toy example
#' 
#' Simulate data for a 3 continuous covariate, 1 binary outcome example where the outcome
#' (\code{Y}) is unrelated to the first covariate (normally-distributed \code{W1}) but is related
#' to the second and third covariates (lognormally-distributed \code{W2} and \code{W3}, which are
#' also related to one another). Note that this function requires the \code{simcausal} package,
#' which may need to be installed from github like so:
#' \code{remotes::install_github('osofr/simcausal', build_vignettes = FALSE)} 
#' 
#' @param n_obs Number of observations to simulate.
#' @param rnd_seed Random seed to be used when simulating.
#' @param verbose Print diagnostic messages? Defaults to FALSE
#' @return A \code{data.frame} with columns \code{ID}, \code{W1}, \code{W2}, \code{W3}, and
#' \code{Y}.
#' @importFrom simcausal simobs
#' @export
#' @examples
#' \dontrun{
#' dat <- sim_toy_data(n_obs = 200, rnd_seed = 620)
#' testobs <- sim_toy_data(n_obs = 1, rnd_seed = 801)
#' }
sim_toy_data <- function(n_obs = 100, rnd_seed = NULL, verbose = FALSE) {
	d <- make_dag_toy_data()
	simobs(d, n = n_obs, rndseed = rnd_seed, verbose = verbose)
}

################################################################################

#' Simulated mixture example with 2 classes.
#' 
#' Simulated mixture example from ElemStatLearn. Two classes. Four + \code{p_noise} covariates.
#' 
#' @param n_obs Number of observations to simulate.
#' @param p_noise Number of N(0,1) random noise covariate columns.
#' @param class_prop Proportion in the second class (labelled 1).
#' @param rnd_seed Random seed to be used when simulating.
#' @return A \code{data.frame}
#' @importFrom MASS mvrnorm
#' @importFrom dplyr bind_rows bind_cols rename_with
#' @importFrom magrittr `%>%`
#' @importFrom tibble remove_rownames
#' @export
#' @examples
#' \dontrun{
#' dat <- sim_elemstatlearn_data(n_obs = 200, rnd_seed = 620)
#' testobs <- sim_elemstatlearn_data(n_obs = 1, rnd_seed = 801)
#' }
sim_elemstatlearn_data <- function(n_obs = 500, p_noise = 6, class_prop = 0.5, rnd_seed = NULL) {

	set.seed(rnd_seed)

	# library(ElemStatLearn)
	# simulated mixture example with 200 instances and two classes. 100 members in each class.
	# data(mixture.example)
	# str(mixture.example)

	Y <- c(rep(0, floor(n_obs*(1-class_prop))), rep(1, ceiling(n_obs*class_prop)))
	# First generate 10 means m_k from a bivariate Gaussian distribution N((1, 0)^T, I) and label this class 0.
	# Similarly, draw 10 more from N((0,1)^T, I) and label class 1.
	means <- list(class0 = mvrnorm(10, mu = c(1,0,1,0), Sigma = diag(4)),
				  class1 = mvrnorm(10, mu = c(0,1,0,1), Sigma = diag(4))) %>%
		lapply(as.data.frame) %>%
		bind_rows()

	# Then for each class generate n/2 observations as follows: 
	#	for each observation, pick an m_k at random with (uniform) probability 1/10, 
	centers <- c(sample(1:10, floor(n_obs*(1-class_prop)), replace = TRUE),
				 sample(11:20, ceiling(n_obs*class_prop), replace = TRUE))
	#	and then generate a N(m_k, I/5), thus leading to a mixture of Gaussian clusters for each class.
	x_mix <- mvrnorm(n_obs, mu = rep(0, 4), Sigma = diag(4)/5) %>%
		as.data.frame() %>%
		rename_with(~ gsub("^V", "Xmix", .x)) %>%
		`+`(means[centers, ])
	x_noise <- mvrnorm(n_obs, mu = rep(0, p_noise), Sigma = diag(p_noise)) %>%
		as.data.frame() %>%
		rename_with(~ gsub("^V", "Xnoise", .x))

	X <- bind_cols(x_mix, x_noise) %>%
		remove_rownames()

	bind_cols(data.frame(ID = seq_along(Y)),
			  X,
			  data.frame(Y = Y))
}

################################################################################

#' Utility function to create indicator variables from a factor
#' 
#' @param col_name Character string; column name in \code{dat} containing a
#' factor variable
#' @param dat \code{data.frame}
#' @return A \code{data.frame} containing indicator column(s)
#' @importFrom stats model.matrix
#' @export
#' @examples
#' factor_to_indicator("Species", iris)
factor_to_indicator <- function(col_name, dat){
    old_na_action <- options('na.action')
    options(na.action = 'na.pass')

    x_factor <- dat[[col_name]]
    x <- gsub("[[:space:]]+", "_", x_factor)
    x_expand <- as.data.frame(model.matrix(~0+x))
    colnames(x_expand) <- paste(col_name, gsub("^x", "", colnames(x_expand)), sep=".")

    if(any(rowSums(is.na(x_expand))>0)) {
        x_miss <- which(rowSums(is.na(x_expand))>0)
        x_expand[x_miss, ] <- 0
        x_miss_col <- rep(0, nrow(x_expand))
        x_miss_col[x_miss] <- 1
        x_expand <- cbind(x_miss_col, x_expand)
        colnames(x_expand)[1] <- paste("miss", col_name, sep = ".")
    }

    # remove baseline level
    baseline_level <- paste(col_name, gsub("[[:space:]]+", "_", levels(x_factor)[1]), sep = ".")
    x_expand <- x_expand[,-which(colnames(x_expand) %in% baseline_level), drop = FALSE]

    options(old_na_action)

    return(x_expand)
}