### Compute post-hoc quantities ###
# Compute densities, means and variances.



#' Approximation to hyperparameter log-posterior
#'
#' @description This function computes log pi(theta|y) for all thetas in the grid.
#'
#' @param optresults Optimization results, a tibble() output by optimize_all_thetas_parallel().
#' @param model_data ccmodeldata object output by model_setup()
#'
#' @return A tibble() in the same format as optresults, with columns added for the log-posterior
#' of theta and sigma.
#'
#' @export
#'
add_log_posterior_values <- function(optresults,model_data) {
  optresults <- dplyr::ungroup(optresults)
  # Log posterior for theta
  logposttheta <- optresults %>%
    purrr::pmap(~log_posterior_theta(unlist(..1),unlist(..4),model_data)) %>%
    as.numeric()
  optresults$theta_logposterior <- logposttheta

  out <- optresults %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sigma = list(exp(-.5 * unlist(.data[["theta"]]))),
                  sigma_logposterior = length(unlist(.data[["sigma"]])) * log(2) - sum(log(unlist(.data[["sigma"]]))) + .data[["theta_logposterior"]])

  attr(out,"thetagrid") <- attributes(optresults)$thetagrid
  out
}

#' Normalize log-posterior via simple numerical integration
#'
#' @description Simple implementation of the trapezoid rule for integration. The key is that
#' this function accepts grid points and LOGGED function values and computes the integral in
#' a numerically stable manner, using matrixStats::logSumExp().
#'
#' For a single hyperparameter, the function works on uniform and non-uniform grids. For
#' a multidimensional hyperparameter, the function only works on a uniform grid (in all dimensions).
#'
#' @param pp logged function values
#' @param tt grid points at which these function values were computed, as a NIGrid object
#'
#' @return A number giving the LOG of the integral of the function used to compute pp.
#' This is the log-normalizing constant.
#'
#' @export
#'
normalize_log_posterior <- function(pp,tt) {
  # tt: grid returned by mvQuad::createNIGrid
  # pp: log posterior evaluated at these points
  ww <- mvQuad::getWeights(tt)
  matrixStats::logSumExp(log(ww) + pp)
}

#' Normalize the log posterior returned by the optimization
#'
#' @description Wrapper around normalize_log_posterior() that takes in
#' the dataframe of optimization results and returns the same dataframe
#' but with the theta_logposterior column normalized.
#'
#' @param optresults Optimization results, a tibble() output by optimize_all_thetas_parallel().
#' @param model_data Model_data needed.
#'
#' @return Data frame of optimization results with the theta_logposterior column normalized,
#' i.e. satisfying logSumExp(theta_logposterior) = 0
#'
#' @export
#'
normalize_optresults_logpost <- function(optresults,model_data) {
  # Get the nodes and weights corresponding to thetas in optresults
  nodesandweights <- cbind(mvQuad::getNodes(model_data$thetagrid),mvQuad::getWeights(model_data$thetagrid))
  K <- ncol(nodesandweights) - 1
  colnames(nodesandweights) <- c(stringr::str_c("theta",1:K),"weights")
  nodesandweights <- as.data.frame(nodesandweights)

  thetaopt <- cbind(purrr::reduce(optresults$theta,rbind),optresults$theta_logposterior)
  colnames(thetaopt) <- c(stringr::str_c("theta",1:K),"theta_logposterior")
  thetaopt <- as.data.frame(thetaopt)

  suppressMessages(thetaoptmerged <- dplyr::left_join(thetaopt,nodesandweights))

  ww <- thetaoptmerged$weights
  pp <- thetaoptmerged$theta_logposterior

  # thetanormconst <- matrixStats::logSumExp(pp + log(ww))
  sp <- matrixStats::logSumExp(pp[ww > 0] + log(ww[ww > 0]))
  sn <- matrixStats::logSumExp(pp[ww < 0] + log(abs(ww)[ww < 0]))
  thetanormconst <- sp - log(1 + exp(sn-sp))
  optresults$theta_logposterior <- optresults$theta_logposterior - thetanormconst

  out <- optresults %>%
    dplyr::rowwise() %>%
    dplyr::mutate(sigma = list(exp(-.5 * unlist(.data[["theta"]]))),
                  sigma_logposterior = length(unlist(.data[["sigma"]])) * log(2) - sum(log(unlist(.data[["sigma"]]))) + .data[["theta_logposterior"]])

  out$weights <- ww

  out
}



#' Marginal means/variances for linear combinations of latent variables
#'
#' @description The compute_marginal_means_and_variances() function lets you specify a sparse matrix of
#' linear combinations of latent variables to compute the marginal means and variances of. The most
#' common use case for this is when you have included an effect as both linear and smooth, and you
#' need the marginal mean and variance of the linear predictor. You need to account for the correlation
#' between the posterior regression coefficient and posterior random effect.
#'
#' Because this is a common use case, this function takes in your model_data and returns a correctly
#' formatted matrix which specifies that you want means/variances for the linear predictors from all
#' terms which appear in the model both as linear and smooth effects.
#'
#' @param model_data ccmodeldata object output by model_setup()
#'
#' @return A sparse matrix with one column per necessary linear combination.
#'
#' @export
#'
make_model_lincombs <- function(model_data) {
  # Check to make sure both linear and smooth terms in the model
  if (length(model_data$model_elements$linear) == 0 | length(model_data$model_elements$smooth) == 0) {
    stop("You should only be looking at linear combinations if there are both linear and smooth terms in your model.")
  }
  # Check to see if the SAME term is included as both linear and smooth.
  if (length(intersect(model_data$model_elements$linear,model_data$model_elements$smooth)) == 0) {
    stop("You have linear and smooth terms in your model, but I don't see any overlap. You should only use this function to create linear combinations for terms in your model that are included as both linear and smooth")
  }

  # Determine the degree of polynomial for each covariate
  polydegrees <- get_polynomial_degree(model_data$model_elements$linear_formula)

  # Get the indices for all terms
  indices <- get_indices(model_data,removezeros = TRUE)
  indices_zeroes <- get_indices(model_data,removezeros = FALSE)
  # Wd <- model_data$Nd + length(indices$smooth) + length(indices$linear)
  Wd <- model_data$Wd
  smooth_terms <- unique(names(indices$smooth))
  linear_terms <- unique(names(indices$linear))
  # The terms we use are the terms that appear both in smooth and linear
  terms_to_use <- intersect(smooth_terms,linear_terms)

  # Create the linear combinations
  lincomblist <- list()
  outmatlist <- list()
  for (nm in terms_to_use) {
    degree <- polydegrees[names(polydegrees) == nm]
    idx <- indices$smooth[names(indices$smooth) == nm]
    linear_idx <- indices$linear[names(indices$linear) == nm]
    u <- indices$covvalues[[nm]]

    for (j in 1:length(u)) {
      betavec <- u[j]^(1:degree)
      ll <- sparseVector(
        x = c(1,betavec),
        i = c(idx[j],linear_idx),
        length = Wd
      )
      lincomblist <- c(lincomblist,ll)
    }
    outmat <- lincomblist %>%
      purrr::map(~as(.,"sparseMatrix")) %>%
      purrr::reduce(cbind)

    # Correct for removed zeros
    if (!is.null(model_data$vectorofcolumnstoremove)) {
      if (length(model_data$vectorofcolumnstoremove) > 0 & !all(model_data$vectorofcolumnstoremove == 0)) {
        constr <- model_data$control$linear_constraints[[nm]]
        wcol <- which(constr$u == constr$whichzero[1])
        nr <- nrow(outmat)
        bindvec <- as(sparseVector(x = constr$u[wcol]^(1:degree),i = linear_idx,length = Wd),"sparseMatrix")
        if (wcol == 1) {
          outmat <- cbind(bindvec,outmat)
        } else if (wcol == ncol(outmat) + 1) {
          outmat <- cbind(outmat,bindvec)
        } else {
          outmat <- cbind(outmat[ ,1:(wcol-1)],bindvec,outmat[ ,wcol:ncol(outmat)])
        }
      }
    }

    outmatlist <- c(outmatlist,outmat)
    lincomblist <- list()
  }

  # Note: the following is not the fastest way to do this. See stackoverflow:
  # https://stackoverflow.com/questions/8843700/creating-sparse-matrix-from-a-list-of-sparse-vectors#8844057
  # It's about 2x - 3x faster in benchmarking; not worth introducing new code.
  outmatlist %>% purrr::reduce(cbind)
}


#' Make a matrix of linear constraints
#'
#' @description Take the linear constraints from the control argument of a ccmodeldata object
#' and create a sparse matrix of linear constraints.
#'
#' @param model_data ccmodeldata object created by model_setup. Must have control$linear_constraints
#'
#' @return A sparse matrix where each column represents one linear constraint. If the return value is A
#' then AW = 0.
#'
#' @export
#'
make_linear_constraints <- function(model_data) {
  if (is.null(model_data$control$linear_constraints)) stop("No linear constraints provided in model_data$control")
  if (length(model_data$control$linear_constraints) == 0) stop("No linear constraints provided in model_data$control")

  # Pull the constraints from the model data, retaining only those that
  # have more than one value. Single-value constraints are automatically done in the precision matrix
  constr <- model_data$control$linear_constraints %>% purrr::map("whichzero") %>%
    purrr::keep(~length(.x) > 1)
  if (length(constr) == 0) return(0) # No additional constraints beyond what was manually set to zero.
  # Get the indices of all unique covariate values that haven't already been set to zero:
  idx <- get_indices(model_data)
  # For the constraints that appear in covvalues (because they haven't already been
  # manually removed), figure out which element of the covvalues it is, then take
  # that element from the idx$smooth for that covariate (!)
  get_constraint_index <- function(nm) {
    cvec <- constr[[nm]]
    ivec <- idx$covvalues[[nm]]
    idxvec <- idx$smooth[names(idx$smooth) == nm]
    cvec <- cvec[cvec %in% ivec]
    out <- numeric(length(cvec))
    for (j in 1:length(cvec)) out[j] <- idxvec[which(cvec[j] == ivec)]
    unname(out)
  }
  # Create a sparse matrix of constraints
  names(constr) %>%
    purrr::map(get_constraint_index) %>%
    purrr::map(~Matrix::sparseVector(x = 1,i = .x,length = model_data$Wd)) %>%
    purrr::map(as,Class = "sparseMatrix") %>%
    purrr::reduce(cbind) %>%
    as(Class = "dgTMatrix")

}


#' Compute marginal means and variances
#'
#' @description Given optimization results with log posteriors computed, compute the
#' marginal means and variances. These are the final outputs of the whole estimation
#' procedure.
#'
#' Linear constraints are corrected for at this point. Marginal variances for linear
#' combinations of latent variables are also available.
#'
#' @param i Either a 1) vector giving the indices of the latent variables for which you want marginal
#' means and variances or 2) an object of class ccindex output by get_indices() which prescribes which
#' terms you want means/variances for.
#' @param model_results Output of optimization; can be before or after you call add_log_posterior_values().
#' @param model_data ccmodeldata object output by model_setup()
#' @param constrA Either a sparse matrix whose columns contain linear constraints under which you would
#' like to compute means/variances, or NULL. If NULL, any linear constraints will be pulled from model_data.
#' @param lincomb Either a sparse matrix whose columns contain linear combinations of the latent variables
#' whose means/variances you would like to compute, or an object of class cclincomb output by make_model_lincombs().
#' If NULL, will be computed automatically. Set lincomb = FALSE in order to prevent this.
#'
#'
#' @export
#'

# i <- index11
# model_results <- opt_11
# model_data <- model_data11

compute_marginal_means_and_variances <- function(i=NULL,model_results,model_data,constrA = NULL,lincomb = NULL) {
  if (nrow(model_results) > 1) {
    # Add log posterior values for theta if not present
    if (!("theta_logposterior" %in% names(model_results))) {
      model_results <- add_log_posterior_values(model_results,model_data)
    }
    # Normalize
    thetanormconst <- normalize_log_posterior(model_results$theta_logposterior,attributes(model_results)$thetagrid)
    model_results$theta_logposterior <- model_results$theta_logposterior - thetanormconst
    # Get the integration weights
    intweights <- mvQuad::getWeights(attributes(model_results)$thetagrid)[ ,1]
  }
  # Compute the precision matrices for each theta
  precision_matrices <- model_results %>%
    purrr::pmap(~list(Q = Q_matrix(theta = ..1,model_data = model_data),
                      theta = ..1)
    )
  # Compute the hessians for each theta
  hessians <- list()
  myhes <- parallel::mclapply(model_results$solution, hessian_log_likelihood,model_data = model_data)
  for (j in 1:length(model_results$theta)) {
    hessians[[j]] <- list(C=myhes[[j]],theta=model_results$theta[j])
  }
  # If linear combinations required, set up the relevant functions
  if (!is.null(lincomb)) {
    compute_var_one_lincomb <- function(a,Q) {
      # a <- cbind(a)
      ZZ <- solve(Q,a)
      as.numeric(crossprod(a,ZZ))
    }
    compute_var_all_lincombs <- function(A,Q) {
      # Coerce to list of sparse vectors
      AA <- list()
      result <- c()
      for (j in 1:ncol(lincomb)){
        AA[[j]] <- as(lincomb[,j],"sparseVector")
        result[j] <- compute_var_one_lincomb(AA[[j]],Q)
      }
      result
    }
    compute_one_lincomb_correction <- function(a,WW,VV) {
      as.numeric(crossprod(crossprod(VV,a),crossprod(WW,a)))
    }
    compute_all_lincomb_correction <- function(A,WW,VV) {
      AA <- list()
      result <- c()
      for (j in 1:ncol(lincomb)){
        AA[[j]] <- as(lincomb[,j],"sparseVector")
        result[j] <- compute_one_lincomb_correction(AA[[j]],WW,VV)
      }
      result
    }
  }
  # If no linear constraints, compute the marginal means and variances as normal
  if (is.null(constrA)) {
    margmeans <- model_results %>%
      purrr::pmap(~..4) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)
    # Marginal variances: add the precision and the hessian and get diagOfInv
    if (is.null(i)) i <- 1:model_data$Wd
    margvars <- purrr::map2(precision_matrices,hessians,~.x[["Q"]] + .y[["C"]]) %>%
      purrr::map(~diag(solve(.x))[i]) %>%
      purrr::reduce(rbind)
    # If there are linear combinations, compute their variances separately from diagOfInv
    if (!is.null(lincomb)) {
      # lincomb is a column matrix. Change to list and map over the columns
      lincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]])) %>%
        purrr::map("lincombvar") %>%
        purrr::reduce(rbind)
    }
  }
  else {
    # If there are linear constraints, compute the corrected mean and variance
    # First compute the uncorrected mean
    uncorrectedmean <- purrr::pmap(model_results,~list(theta = ..1,mode = ..4))
    # Get the precision matrix of the GMRF- Q + C
    QpC <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]]))
    # Compute the correction term
    WW <- purrr::map(QpC,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
    YY <- purrr::map2(WW,uncorrectedmean,~list(
      YY = solve(t(constrA) %*% .x[["WW"]],t(constrA) %*% .y[["mode"]]),
      theta = .y[["theta"]]
    ))
    correction_mean <- purrr::map2(WW,YY,~list(correction = .x[["WW"]] %*% .y[["YY"]],theta = .x[["theta"]]))

    # Now correct
    margmeans <- purrr::map2(uncorrectedmean,correction_mean,
                             ~.x[["mode"]] - .y[["correction"]]) %>%
      purrr::map(t) %>%
      purrr::reduce(rbind)
    # Now compute the variances
    # Add the corrected mean to the model_results
    model_results$corrected_mean <- vector(mode = "list",length = nrow(model_results))
    for (k in 1:nrow(model_results)) model_results$corrected_mean[[k]] <- margmeans[k, ]
    # Re-compute the hessians
    corrected_hessians <- list()
    myhes2 <- parallel::mclapply(model_results$corrected_mean, hessian_log_likelihood,model_data = model_data)
    for (j in 1:length(model_results$theta)) {
      corrected_hessians[[j]] <- list(C=myhes2[[j]],theta=model_results$theta[j])
    }
    # Get the corrected precision matrix of the GMRF- Q + C_correct
    QpC_corrected <- purrr::map2(precision_matrices,corrected_hessians,~list(QpC = as((as.matrix(.x[["Q"]] + .y[["C"]])),"sparseMatrix"),theta = .x[["theta"]]))
    # uncorrectedvariances <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = NULL,i = i))
    margvars <- purrr::map(QpC_corrected,~diagOfInv(x = .x[["QpC"]],constrA = constrA,i = i)) %>%
      purrr::reduce(rbind)


    if (!is.matrix(margvars)) margvars <- matrix(margvars,nrow = 1)
    # If we require marginal variances for linear combinations, compute them separately
    if (!is.null(lincomb)) {
      uncorrectedlincombvars <- purrr::map2(precision_matrices,hessians,~list(QpC = .x[["Q"]] + .y[["C"]],theta = .x[["theta"]])) %>%
        purrr::map(~list(lincombvar = compute_var_all_lincombs(lincomb,.x[["QpC"]]),theta = .x[["theta"]]))

      # Compute the corrections
      WW <- purrr::map(QpC_corrected,~list(WW = solve(.x[["QpC"]],constrA),theta = .x[["theta"]]))
      VV <- purrr::map(WW,~list(VV = solve(t(.x[["WW"]]) %*% constrA,t(.x[["WW"]])),theta = .x[["theta"]])) %>%
        purrr::map(~list(VV = t(.x[["VV"]]),theta = .x[["theta"]]))

      lincombvarcorrections <- list()
      for (jj in 1:length(WW)) {
        lincombvarcorrections[[jj]] <- list(lincombvarcorrection = compute_all_lincomb_correction(lincomb,WW[[jj]]$WW,VV[[jj]]$VV),theta = WW[[jj]]$theta)
      }

      lincombvars <- purrr::map2(uncorrectedlincombvars,lincombvarcorrections,
                                 ~list(lincombvar = .x[["lincombvar"]] - .y[["lincombvarcorrection"]],
                                       theta = .x[["theta"]])) %>%
        purrr::map("lincombvar") %>%
        purrr::reduce(rbind)
    }
  }
  if (nrow(margmeans) == 1) {
    finalmeans <- as.numeric(margmeans)[i]
    finalvars <- as.numeric(margvars)
    finallincombvars <- NULL
    if (!is.null(lincomb))
      finallincombvars <- as.numeric(lincombvars)
  }
  else {
    postvals <- exp(model_results$theta_logposterior + log(intweights))
    finalmeans <- sweep(margmeans,1,postvals,"*") %>% apply(2,sum)
    finalvars <- sweep(margvars,1,postvals,"*") %>% apply(2,sum)
    finallincombvars <- NULL
    if (!is.null(lincomb)) finallincombvars <- sweep(lincombvars,1,postvals,"*") %>% apply(2,sum)
    finalmeans <- finalmeans[i]
  }
  list(mean = finalmeans,
       variance = finalvars,
       lincombvars = finallincombvars)
}


marginal_hyperparameter_posterior <- function(j,optresults,quantiles = c(2.5,97.5)/100, model_data){
  thetagridfull <- model_data$thetagrid
  S <- thetagridfull$dim
  # If it's already one-dimensional, don't need to do anything new, but do compute quantiles
  if (S == 1) {
    outmat <- dplyr::tibble(
      theta = purrr::reduce(optresults$theta,rbind),
      thetalogmargpost = optresults$theta_logposterior,
      sigma = purrr::reduce(optresults$sigma,rbind),
      sigmalogmargpost = optresults$sigma_logposterior
    )

    thetacumsum <- cumsum(mvQuad::getWeights(thetagridfull) * exp(outmat$thetalogmargpost))
    thetaquantiles <- purrr::map(quantiles,~outmat$theta[which(thetacumsum == min(thetacumsum[thetacumsum > .x]))]) %>% as.numeric()
    sigmaquantiles <- rev(exp(-.5 * thetaquantiles))

    return(
      list(
        margpost = outmat,
        quantiles = dplyr::tibble(whichmarginal = rep(1,length(quantiles)),q = quantiles,theta = thetaquantiles,sigma = sigmaquantiles)
      )
    )
  }
  # Get the reduced grid
  thetagridreduced <- mvQuad::createNIGrid(
    dim = thetagridfull$dim - 1,
    type = thetagridfull$type[-j],
    level = as.numeric(thetagridfull$level[ ,-j]),
    ndConstruction = thetagridfull$ndConstruction,
    level.trans = thetagridfull$level.trans
  )
  mvQuad::rescale(thetagridreduced,domain = thetagridfull$features$domain[-j, ])

  # Get a 1-d grid, for computing quantiles at the end.
  thetagrid1d <- mvQuad::createNIGrid(
    dim = 1,
    type = thetagridfull$type[j],
    level = as.numeric(thetagridfull$level[ ,j]),
    ndConstruction = thetagridfull$ndConstruction,
    level.trans = thetagridfull$level.trans
  )
  mvQuad::rescale(thetagrid1d,domain = thetagridfull$features$domain[j, ])

  # Return the marginal posterior and its evaluation points
  # In the optimization results, we have a matrix of theta values which matches the full grid,
  # and the log posterior evaluated at these values.
  nodesfull <- purrr::reduce(optresults$theta,rbind)
  nodesfull <- cbind(nodesfull,optresults$theta_logposterior) # Theta logposterior is now the last column
  nodesfull <- nodesfull[order(nodesfull[ ,j]), ]
  colnames(nodesfull) <- c(paste0("theta",1:S),"thetalogpost")
  nodesfull <- as.data.frame(nodesfull)

  # Now we have a matrix of thetavalues, nodes, and function values-- add on the weights
  nodesmulti <- mvQuad::getNodes(thetagridreduced)
  nodesmulti <- cbind(nodesmulti,mvQuad::getWeights(thetagridreduced))
  colnames(nodesmulti) <- c(paste0("theta",(1:S)[-j]),"weights")
  nodesmulti <- as.data.frame(nodesmulti)

  suppressMessages({# It prints what it's joining by, which is all columns, and I don't want to see this printed
    thetamargposts <- dplyr::left_join(nodesfull,nodesmulti) %>%
      dplyr::group_by(.data[[stringr::str_c("theta",j)]]) %>%
      dplyr::summarize(thetalogmargpost = matrixStats::logSumExp(.data[["thetalogpost"]] + log(.data[["weights"]])))
  })

  thetamargposts$whichmarginal <- rep(j,nrow(thetamargposts))
  # Now add on the sigmas
  outmat <- thetamargposts %>%
    dplyr::mutate(sigma = exp(-.5 * .data[[paste0("theta",j)]]),
                  sigmalogmargpost = log(2/.data[["sigma"]]) + .data[["thetalogmargpost"]]
    ) %>%
    dplyr::rename(theta = .data[[paste0("theta",j)]])

  # Mean and sd
  ww <- mvQuad::getWeights(thetagridfull)[ ,1]
  thetamean <- apply(ww * purrr::reduce(optresults$theta,rbind) * exp(optresults$theta_logposterior),2,sum)
  sigmamean <- apply(ww * purrr::reduce(optresults$sigma,rbind) * exp(optresults$theta_logposterior),2,sum)
  thetasd <- sqrt(apply(ww * (purrr::reduce(optresults$theta,rbind) - thetamean)^2 * exp(optresults$theta_logposterior),2,sum))
  sigmasd <- sqrt(apply(ww * (purrr::reduce(optresults$sigma,rbind) - sigmamean)^2 * exp(optresults$theta_logposterior),2,sum))


  # Quantiles
  thetacumsum <- cumsum(mvQuad::getWeights(thetagrid1d) * exp(outmat$thetalogmargpost))
  thetaquantiles <- purrr::map(quantiles,~outmat$theta[min(which(thetacumsum == min(thetacumsum[thetacumsum > .x])))]) %>% purrr::reduce(c)
  sigmaquantiles <- rev(exp(-.5 * thetaquantiles))

  list(
    margpost = outmat,
    margmoments = dplyr::tibble(moment = c("mean","sd"),theta = c(thetamean[j],thetasd[j]),sigma = c(sigmamean[j],sigmasd[j])),
    quantiles = dplyr::tibble(whichmarginal = rep(j,length(quantiles)),q = quantiles,theta = thetaquantiles,sigma = sigmaquantiles)
  )
}


