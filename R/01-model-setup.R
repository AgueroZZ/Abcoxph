### MODEL SETUP ###
# This file contains functions which take in data and
# model specifications and return objects usable
# by the internal likelihood/optimization/summary functions.
abcox_control <- function(...) {
  outargs <- abcox_default_control()
  if (length(list(...)) == 0) return(outargs)
  for (nm in names(list(...))) {
    outargs[[nm]] <- list(...)[[nm]]
  }
  outargs
}


#' Do the AGHQ to find the theta_grid
#'
#' @description Implement AGHQ to optimize for theta_grid
#'
#' @return An NIGrid object of theta_grid
#'
#' @param model_data the model_data created by abcoxph_setup or abcoxph_AGHQ_setup
#' @param k number of points in the grid
#' @param startingvals the starting value of theta in optimization
#'
#'
useAGHQ <- function(model_data,k = 3, startingvals = 0, method = "BFGS", inner_method = "trust"){
  log_posterior_joint <- function(W,theta){
    Q <- Q_matrix(theta,model_data = model_data)
    model_data$theta_logprior(theta) + logprior_W_aghq(W,theta,Q,model_data = model_data) + log_likelihood(W,model_data = model_data)
  }
  grad_log_posterior_W_opt <- function(W,theta,Q = NULL) {
    if (is.null(Q)) {
      Q <- Q_matrix(theta,model_data = model_data)
    }
    - as.numeric(crossprod(Q,W)) + grad_log_likelihood(W,model_data = model_data)
  }
  hessian_log_posterior_W_opt <- function(W,theta = NULL,Q = NULL) {
    if (is.null(theta) & is.null(Q)) stop("One of Q or theta must be provided")
    if (is.null(Q)) Q <- Q_matrix(theta, model_data = model_data)
    A <- as.matrix(-(Q + hessian_log_likelihood(W, model_data = model_data)))
    A <- as(A,'dgCMatrix')
    A
  }
  ff <- list(
    fn = log_posterior_joint,
    gr = grad_log_posterior_W_opt,
    he = function(W,theta) hessian_log_posterior_W_opt(W,theta)
  )
  result <- aghq::marginal_laplace(ff,k,list(W = rep(0,model_data$Wd),theta = startingvals), control = list(method = c(method), inner_method = c(inner_method)))
}






#' Set up the prior distributions that you want to use.
#'
#' @description Setting up your prior distributions to be used.
#'
#' @return A list that contains all the information needed for your prior to be fitted using abcoxph_setup.
#'
#' @param theta_dim How many hyperparameters do you have?
#' @param theta_number For each parameter, how many grids you want to use for the integration?
#' @param theta_grid Your specific grid for theta to be used.
#' @param theta_prior_a A values for the pc.prior, the probability
#' @param theta_prior_u A values for the pc.prior, the lower bound
#' @param beta_logprec A values for the log precision of you fixed effect parameter.
#'
prior_setup <- function(theta_dim = 1, theta_number = 50, theta_grid = c(-2,18), theta_prior_a = 0.5, theta_prior_u = 2, beta_logprec = log(0.001)){
  prior_control <- list()
  prior_control$theta$dim <- theta_dim
  prior_control$theta$num <- theta_number
  prior_control$theta$grid <- theta_grid
  prior_control$theta$prior_a <- theta_prior_a
  prior_control$theta$prior_u <- theta_prior_u
  prior_control$beta_logprec <- beta_logprec
  prior_control
}




#' Set up the prior distributions that you want to use when you use AGHQ setting.
#'
#' @description Setting up your prior distributions to be used, when you plan to use AGHQ in your inference.
#'
#' @return A list that contains all the information needed for your prior to be fitted using abcoxph_AGHQ_setup
#'
#' @param theta_dim How many hyperparameters do you have?
#' @param theta_number For each hyperparameter, how many grids you want to use for the integration?
#' @param startingvals The value where you start the optimization for theta
#' @param theta_prior_a A values for the pc.prior, the probability
#' @param theta_prior_u A values for the pc.prior, the lower bound
#' @param beta_logprec A values for the log precision of you fixed effect parameter.
#'
prior_setup_AGHQ <- function(theta_dim = 1, theta_number = 3, startingvals = 0, theta_prior_a = 0.5, theta_prior_u = 2, beta_logprec = log(0.001)){
  prior_control <- list()
  prior_control$theta$dim <- theta_dim
  prior_control$theta$num <- theta_number
  prior_control$theta$prior_a <- theta_prior_a
  prior_control$theta$prior_u <- theta_prior_u
  prior_control$beta_logprec <- beta_logprec
  prior_control$startingvals <- startingvals
  prior_control
}







#' Set up the Cox Proportional Hazard Model that you want to infer with AGHQ method.
#'
#' @description Setting up a target CoxPH Model and its AGHQ specification, by specifying the formula. The response in the variable should be your
#' observed times, and the regressors's type are well identified with s denote smoothing and id denote between-subject
#' index, other terms are assumed to be linear fixed effect.
#'
#' @return A list that contains all the information needed for your model to be fitted using abcoxph_fit.
#'
#' @param formula A formula specifies your observed times as response variable, regressors type are smoothing, id, or linear fixed effect.
#' @param cens A STRING specifies the name of your censoring indicator in the dataframe.
#' @param data Target dataframe.
#' @param prior_control A list specifies your prior distributions, this must be generated from prior_setup_AGHQ.
#' @param RW2BINS A number specifies the number of bins you want to use for smoothing.
#' @param method The outer method used in AGHQ for the optimization, suggested to use "BFGS" at all the time.
#' @param inner_method The inner method used in AGHQ for the optimization, suggested to use "trust", if the dataset is large, suggested to use "SR1"
#' @param correction Whether to use the Breslow's method for the correction of ties, suggested when the dataset is small and number of ties large.
#'
abcoxph_AGHQ_setup <- function(formula,cens,data,prior_control, RW2BINS = NULL, method = "BFGS", inner_method = "trust", correction = "FALSE"){
  str <- parse_formula(formula)
  model_data <- list(M = 0, p = length(str$linear), n = nrow(data), Nd = nrow(data) - 1)
  model_data$correction <- correction
  if(length(str$linear) > 0){
    newdata <- data[str$linear]
    newdata$times <- data[[str$response]]
  }
  else{
    newdata <- data_frame(times = data[[str$response]])
  }
  newdata$censoring <- data[[cens]]
  newdata$entry <- rep(0,length(newdata$times))
  if(length(str$id) != 0){
    newdata$id <- data[[str$id]]
  }
  if(length(str$smooth) != 0){
    if(is.null(RW2BINS)) RW2BINS <- 50
    else RW2BINS <- RW2BINS
    model_data$RW2BINS <- RW2BINS
    newdata[[str$smooth]] <- data[[str$smooth]]
    newdata[[str$smooth]] <- bin_covariate(newdata[[str$smooth]],bins = RW2BINS,type = "quantile")
  }
  newdata <- arrange_data(newdata)
  newdata$od <- 1:nrow(newdata)
  model_data$times <- newdata$times
  model_data$censoring <- newdata$censoring
  model_data$entry <- newdata$entry
  model_data$od <- newdata$od
  xdata <- newdata[str$linear]
  model_data$X <- as.matrix(xdata)
  model_data$diffmat <- create_diff_matrix(model_data$n)
  model_data$lambdainv <- Matrix::diag(model_data$n-1)
  model_data$Xd <- model_data$diffmat %*% model_data$X
  if(length(str$id) != 0){
    Blist <- list()
    Blist$id <- create_blist_element(newdata$id)
    iidnum <- Blist %>% map("B") %>% map(ncol) %>% reduce(sum)
    model_data$B <- Blist
    model_data$M <- model_data$M + iidnum
    model_data$B$id$Bd <- model_data$diffmat %*% model_data$B$id$B
    model_data$modelspec <- model_data$B %>%
      purrr::map("model") %>%
      purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
      purrr::reduce(bind_rows)
  }
  if(length(str$smooth) != 0){
    Alist <- list()
    Alist[[str$smooth]] <- create_alist_element(newdata[[str$smooth]])
    rw2num <- Alist %>% map("A") %>% map(ncol) %>% reduce(sum) - 1
    model_data$A <- Alist
    model_data$M <- model_data$M + rw2num
    model_data$A[[str$smooth]][["Ad"]] <- model_data$diffmat %*% model_data$A[[str$smooth]][["A"]]
    model_data$modelspec <- model_data$A %>%
      purrr::map("model") %>%
      purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
      purrr::reduce(bind_rows)
    model_data$vectorofcolumnstoremove <- round(RW2BINS/2)
  }
  model_data$Ne <- model_data$n
  model_data$Wd <- model_data$M + model_data$p + model_data$Nd
  model_data$Wdf <- model_data$M + model_data$p + model_data$Ne
  model_data$beta_logprec <- prior_control$beta_logprec
  model_data$theta_logprior <- function(theta,prior_alpha = prior_control$theta$prior_a, prior_u = prior_control$theta$prior_u) {
    # In this model, theta is the LOG PRECISION of the rw2 smoothing variance
    # Implement the PC prior directly.
    # P(sigma > u) = alpha.
    # See inla.doc("pc.prec")
    lambda <- -log(prior_alpha)/prior_u
    log(lambda/2) - lambda * exp(-theta/2) - theta/2
  }
  model_data$thetaAGHQ <- useAGHQ(model_data,k = prior_control$theta$num,startingvals = prior_control$startingvals, method = method, inner_method = inner_method)
  model_data$thetagrid <- model_data$thetaAGHQ$normalized_posterior$grid
  model_data
}





#' Set up the Cox Proportional Hazard Model that you want to infer.
#'
#' @description Setting up a target CoxPH Model, by specifying the formula. The response in the variable should be your
#' observed times, and the regressors's type are well identified with s denote smoothing and id denote between-subject
#' index, other terms are assumed to be linear fixed effect.

#'
#' @return A list that contains all the information needed for your model to be fitted using abcoxph_fit.
#'
#' @param formula A formula specifies your observed times as response variable, regressors type are smoothing, id, or linear fixed effect.
#' @param cens A STRING specifies the name of your censoring indicator in the dataframe.
#' @param data Target dataframe.
#' @param prior_control A list specifies your prior distributions.
#' @param RW2BINS A number specifies the number of bins you want to use for smoothing.
#'
abcoxph_setup <- function(formula,cens,data,prior_control = NULL, RW2BINS = NULL){
  str <- parse_formula(formula)
  model_data <- list(M = 0, p = length(str$linear), n = nrow(data), Nd = nrow(data) - 1)
  if(length(str$linear) > 0){
    newdata <- data[str$linear]
    newdata$times <- data[[str$response]]
  }
  else{
    newdata <- data_frame(times = data[[str$response]])
  }
  newdata$censoring <- data[[cens]]
  newdata$entry <- rep(0,length(newdata$times))
  if(length(str$id) != 0){
    newdata$id <- data[[str$id]]
  }
  if(length(str$smooth) != 0){
    if(is.null(RW2BINS)) RW2BINS <- 50
    else RW2BINS <- RW2BINS
    model_data$RW2BINS <- RW2BINS
    newdata[[str$smooth]] <- data[[str$smooth]]
    newdata[[str$smooth]] <- bin_covariate(newdata[[str$smooth]],bins = RW2BINS,type = "quantile")
  }
  newdata <- arrange_data(newdata)
  newdata$od <- 1:nrow(newdata)
  model_data$times <- newdata$times
  model_data$censoring <- newdata$censoring
  model_data$entry <- newdata$entry
  model_data$od <- newdata$od
  xdata <- newdata[str$linear]
  model_data$X <- as.matrix(xdata)
  model_data$diffmat <- create_diff_matrix(model_data$n)
  model_data$lambdainv <- Matrix::diag(model_data$n-1)
  model_data$Xd <- model_data$diffmat %*% model_data$X
  if(length(str$id) != 0){
    Blist <- list()
    Blist$id <- create_blist_element(newdata$id)
    iidnum <- Blist %>% map("B") %>% map(ncol) %>% reduce(sum)
    model_data$B <- Blist
    model_data$M <- model_data$M + iidnum
    model_data$B$id$Bd <- model_data$diffmat %*% model_data$B$id$B
    model_data$modelspec <- model_data$B %>%
      purrr::map("model") %>%
      purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
      purrr::reduce(bind_rows)
  }
  if(length(str$smooth) != 0){
    Alist <- list()
    Alist[[str$smooth]] <- create_alist_element(newdata[[str$smooth]])
    rw2num <- Alist %>% map("A") %>% map(ncol) %>% reduce(sum) - 1
    model_data$A <- Alist
    model_data$M <- model_data$M + rw2num
    model_data$A[[str$smooth]][["Ad"]] <- model_data$diffmat %*% model_data$A[[str$smooth]][["A"]]
    model_data$modelspec <- model_data$A %>%
      purrr::map("model") %>%
      purrr::map2(.,names(.),~tibble(covariate = .y,model = .x)) %>%
      purrr::reduce(bind_rows)
    model_data$vectorofcolumnstoremove <- round(RW2BINS/2)
  }
  model_data$Ne <- model_data$n
  model_data$Wd <- model_data$M + model_data$p + model_data$Nd
  model_data$Wdf <- model_data$M + model_data$p + model_data$Ne

  #### Prior Setting
  if(length(prior_control) == 0){
    model_data$thetagrid <- mvQuad::createNIGrid(dim = 1,type = "GLe",level = 60)
    mvQuad::rescale(model_data$thetagrid,domain = c(-1,8))
    model_data$beta_logprec <- log(.001)
    model_data$theta_logprior <- function(theta,prior_alpha = c(.5),prior_u = c(2)) {
      # In this model, theta is the LOG PRECISION of the rw2 smoothing variance
      # Implement the PC prior directly.
      # P(sigma > u) = alpha.
      # See inla.doc("pc.prec")
      lambda <- -log(prior_alpha)/prior_u
      log(lambda/2) - lambda * exp(-theta/2) - theta/2
    }
  }
  else{
    model_data$thetagrid <- mvQuad::createNIGrid(dim = prior_control$theta$dim,type = "GLe",level = prior_control$theta$num)
    mvQuad::rescale(model_data$thetagrid,domain = prior_control$theta$grid)
    model_data$beta_logprec <- prior_control$beta_logprec
    model_data$theta_logprior <- function(theta,prior_alpha = prior_control$theta$prior_a, prior_u = prior_control$theta$prior_u) {
      # In this model, theta is the LOG PRECISION of the rw2 smoothing variance
      # Implement the PC prior directly.
      # P(sigma > u) = alpha.
      # See inla.doc("pc.prec")
      lambda <- -log(prior_alpha)/prior_u
      log(lambda/2) - lambda * exp(-theta/2) - theta/2
    }
  }
  model_data
}


abcox_default_control <- function(){
  res <- list(
    prec = 1e-08,
    stop.trust.radius = 1e-10,
    report.freq = 10,
    report.level = 4,
    start.trust.radius = 100,
    contract.threshold = .25,
    contract.factor = .5,
    expand.factor = 3,
    trust.iter = 2000,
    maxit = 3000,
    preconditioner = 0
  )
  res
}

#' Do the Approximate Bayesian Inference for the specified Cox Proportional Hazard Model.
#'
#' @description Using the method described in the paper, to perform the approximate Bayesian inference for the specified CoxPH model,
#' under you specified controls.
#'
#' @return A list that contains marginal_latent, marginal_hyper and the original model_data from your model.
#'
#' @param abcox_model A result from abcoxph_setup
#' @param control A list from abcox_default_control
#' @param PARALLEL_EXECUTION Whether you want the computations to be done in parallel.
#'
abcox_fit <- function(abcox_model,control = abcox_default_control(), PARALLEL_EXECUTION = T, method = ""){
  sim1opt <- optimize_all_thetas_parallel(
    theta = abcox_model$thetagrid,
    model_data = abcox_model,
    optcontrol = control,
    doparallel = PARALLEL_EXECUTION
  )
  optresults_withlogpostpre <- add_log_posterior_values(sim1opt,model_data = abcox_model)
  optresults_withlogpost <- normalize_optresults_logpost(optresults_withlogpostpre, model_data = abcox_model)
  margpost1 <- marginal_hyperparameter_posterior(1,optresults_withlogpost, model_data = abcox_model)
  margmeans_and_vars <- compute_marginal_means_and_variances(
    i = (abcox_model$Nd + 1):(abcox_model$Wd),
    model_results = optresults_withlogpostpre,
    model_data = abcox_model,
    lincomb = NULL,
    constrA = NULL
  )
  if (length(abcox_model$RW2BINS) != 0) {
    #smooth
    margmeans_rw2 <- margmeans_and_vars$mean[1:(abcox_model$RW2BINS - 1)]
    margvars_rw2 <- margmeans_and_vars$variance[1:(abcox_model$RW2BINS - 1)]
    margsd_rw2 <- sqrt(margvars_rw2)
    #iid
    margmeans_iid <- NULL
    margvars_iid <- NULL
    margsd_iid <- NULL
    if ((abcox_model$M) >= abcox_model$RW2BINS){
      margmeans_iid <- margmeans_and_vars$mean[abcox_model$RW2BINS:(abcox_model$M)]
      margvars_iid <- margmeans_and_vars$variance[abcox_model$RW2BINS:(abcox_model$M)]
      margsd_iid <- sqrt(margvars_iid)
    }
    #linear
    margmeans_linear <- NULL
    margsd_linear <- NULL
    if ((abcox_model$p) > 0){
    margmeans_linear <- margmeans_and_vars$mean[(abcox_model$M + 1):length(margmeans_and_vars$mean)]
    margsd_linear <- sqrt(margmeans_and_vars$variance[(abcox_model$M + 1):length(margmeans_and_vars$mean)])
    }
    vv <- abcox_model$vectorofcolumnstoremove
    margmeans_rw2 <- c(
      margmeans_rw2[1:(vv - 1)],
      0,
      margmeans_rw2[vv:length(margmeans_rw2)]
    )
    margsd_rw2 <- c(
      margsd_rw2[1:(vv - 1)],
      0,
      margsd_rw2[vv:length(margsd_rw2)]
    )
  }else if (abcox_model$M != 0) {
    margmeans_iid <- margmeans_and_vars$mean[1:(abcox_model$M)]
    margvars_iid <- margmeans_and_vars$variance[1:(abcox_model$M)]
    margsd_iid <- sqrt(margvars_iid)
    margmeans_rw2 <- NULL
    margsd_rw2 <- NULL
    margmeans_linear <- margmeans_and_vars$mean[(abcox_model$M + 1):length(margmeans_and_vars$mean)]
    margsd_linear <- sqrt(margmeans_and_vars$variance[(abcox_model$M + 1):length(margmeans_and_vars$mean)])
  }
  else{
    margmeans_iid <- NULL
    margsd_iid <- NULL
    margmeans_rw2 <- NULL
    margsd_rw2 <- NULL
    margmeans_linear <- margmeans_and_vars$mean
    margsd_linear <- sqrt(margmeans_and_vars$variance)
  }
  marginal_latent <- list(marginal_mean_smoothing = margmeans_rw2, marginal_sd_smoothing = margsd_rw2, marginal_mean_iid = margmeans_iid, marginal_sd_iid = margsd_iid,marginal_mean_linear =  margmeans_linear, marginal_sd_linear = margsd_linear)
  marginal_hyper <- margpost1
  list(marginal_latent = marginal_latent,marginal_hyper = marginal_hyper, theta_logprior = abcox_model$theta_logprior, model_data = abcox_model)
}

#' Plot the posterior for theta.
#'
#' @description Plot the posterior for theta i.e.log-precision, which is attained using laplace approximation.
#'
#' @return A Plot for theta's posterior.
#'
#' @param abcox_fit A result from abcox_fit
#'
plot_theta <- function(abcox_fit){
  priorfunc <- function(x) exp(abcox_fit$theta_logprior(x))
  margpost1 <- abcox_fit$marginal_hyper
  thetapostplot1 <- margpost1$margpost %>%
    mutate(theta_post = exp(thetalogmargpost)) %>%
    ggplot2::ggplot(ggplot2::aes(x = theta)) +
    ggplot2::theme_classic() +
    ggplot2::geom_line(ggplot2::aes(y = theta_post),colour = "black", linetype = "solid", size = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = priorfunc(theta)),colour = "black", linetype = 'dashed', size = 0.5) +
    ggplot2::labs(x = "",y = "") +
    ggplot2::theme_classic(base_size = 28)
  thetapostplot1
}

#' Plot the posterior for sigma
#'
#' @description Plot the posterior for sigma i.e.the standard deviation, which is attained using laplace approximation with
#' a transformation.
#'
#' @return A Plot for sigma's posterior.
#'
#' @param abcox_fit A result from abcox_fit
#'
plot_sigma <- function(abcox_fit){
  priorfuncsigma <- function(x) (2/x) * exp(abcox_fit$theta_logprior(-2*log(x)))
  margpost1 <- abcox_fit$marginal_hyper
  sigmapostplot1 <- margpost1$margpost %>%
    mutate(sigma_post = exp(sigmalogmargpost)) %>%
    ggplot2::ggplot(ggplot2::aes(x = sigma)) +
    ggplot2::theme_classic() +
    ggplot2::geom_line(ggplot2::aes(y = sigma_post),colour = "black",linetype = "solid",size = 0.5) +
    ggplot2::labs(x = "",y = "") +
    ggplot2::geom_line(ggplot2::aes(y = priorfuncsigma(sigma)),colour = 'black',linetype = 'dashed',size = 0.5) +
    ggplot2::theme_classic(base_size = 28)
  sigmapostplot1
}

#' Plot the inferred smoothing function.
#'
#' @description Plot the posterior mean of the smoothing function, and its corresponding credible interval.
#'
#' @return A Plot for the smoothing result.
#'
#' @param abcox_fit A result from abcox_fit
#'
#'
plot_smoothing <- function(abcox_fit){
  model_data <- abcox_fit$model_data
  margmeanall <- abcox_fit$marginal_latent$marginal_mean_smoothing
  margsd <- abcox_fit$marginal_latent$marginal_sd_smoothing
  simplot <- dplyr::tibble(
    x = sort(unique(model_data$A[[1]]$u)),
    mymean = margmeanall,
    mymeanlower = mymean - 2*margsd,
    mymeanupper = mymean + 2*margsd
  ) %>%
    ggplot2::ggplot(ggplot2::aes(x = x)) +
    ggplot2::theme_light() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = exp(mymeanlower),ymax = exp(mymeanupper)),fill = "lightgrey",alpha = .5) +
    ggplot2::geom_line(ggplot2::aes(y = exp(mymeanupper)),colour = "black",linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(y = exp(mymeanlower)),colour = "black",linetype = "dotted") +
    ggplot2::geom_line(ggplot2::aes(y = exp(mymean)),colour = 'black',linetype = 'dotdash') +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme_classic(base_size = 28)
  simplot
}



