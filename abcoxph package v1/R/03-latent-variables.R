### LATENT FIELD ###

# Functions to implement calculations relating to the latent variables, priors and posteriors.
# Most of these are not exported.


Q_matrix_linear <- function(theta,model_data,tau = exp(7), debug = FALSE) {
  # Compute the Q matrix for a model with only linear terms
  # Arguments:
  #   theta: scalar, value of LOG PRECISION log(1/variance) of prior distribution on beta
  #   model_data: usual model data object. Contains the differenced design matrix Xd
  #   tau: INLA's linear predictor noise term
  # Returns:
  #   A large, sparse matrix representing Q from the paper.

  # UPDATE: theta is obsolete. will pull this from model data
  if ("beta_logprec" %in% names(model_data)) theta <- model_data$beta_logprec

  # Construct the matrix
  Sbinv <- Matrix::Diagonal(ncol(model_data$Xd),exp(theta))
  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% model_data$Xd),
    cbind(-tau*t(model_data$Xd) %*% model_data$lambdainv,Sbinv + tau*Matrix::crossprod(Matrix::crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd))
  )
}


Q_matrix_iid_one_component <- function(theta,model_data,covariate) {
  # theta: log(precision); log of 1/(iid random effect variance)
  # covariate: character, name of covariate as appears in model_data$B
  # This function creates the component of the Q matrix corresponding to a single covariate.
  # model_data contains an element "B" which is a list; names(B) is the vector of names of
  # covariates to be modelled as iid random effect. B$covariate is itself a list containing Bmat, the
  # random effects design matrix, and u, the vector containing the UNIQUE values of the covariate. The columns
  # of B$covariate$Bmat should already be ordered according to order(unique(u)); u can thus be either
  # sorted or not. It will be sorted inside this function
  # u should ALREADY be unique and sorted. This just makes this explicit in the code.
  theta2 <- theta[length(theta)]
  u <- sort(unique(model_data$B[[covariate]]$u))
  ul <- length(u)
  AA <- Matrix::Diagonal(n = ul,1)
  exp(theta2) * AA
}

Q_matrix_iid <- function(theta,model_data,tau = exp(7)) {
  # Figure out how many rw2 components there are
  if (is.null(model_data$B)) stop("no iid components in model")

  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)

  howmanyiid <- length(whichiid)


  Suinv <- purrr::map2(whichiid,1:howmanyiid,
                       ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)

  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% Bd),
    cbind(-tau*t(Bd)%*% model_data$lambdainv,Suinv + tau * Matrix::crossprod(Bd,Matrix::crossprod(model_data$lambdainv,Bd)))
  )
}





Q_matrix_rw2_one_component <- function(theta,model_data,covariate) {
  theta1 <- theta[1]
  u <- sort(unique(model_data$A[[covariate]]$u))

  ul <- length(u)
  du <- diff(u)

  H <- Matrix::bandSparse(n = ul,
                          diagonals = list(
                            c(1/du[-(ul-1)],0),
                            c(0,-(1/du[-(ul-1)] + 1/du[-1]),0),
                            c(0,1/du[-1])
                          ),
                          k = c(-1,0,1))

  AA <- Matrix::Diagonal(n = ul,x = c(2/du[1],2/(du[-(ul-1)] + du[-1]),2/du[(ul-1)]))

  exp(theta1) * Matrix::forceSymmetric(Matrix::crossprod(H,Matrix::crossprod(AA,H)))
}



Q_matrix_rw2 <- function(theta,model_data,tau = exp(7), buffer = 1/exp(7)) {
  # Figure out how many rw2 components there are
  if (is.null(model_data$A)) stop("no rw2 components in model")

  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)

  howmanyrw2 <- length(whichrw2)


  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% map("Ad") %>% purrr::reduce(cbind)

  # APPLY ZERO CONSTRAINTS
  # Directly apply constraint that U_t = 0 by deleting the t^th column of A,
  # and the t^th row and column of Suinv

  if (!(0 %in% model_data$vectorofcolumnstoremove)) {
    Ad <- Ad[ ,-model_data$vectorofcolumnstoremove]
    Suinv <- Suinv[-model_data$vectorofcolumnstoremove,-model_data$vectorofcolumnstoremove]
  }

  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(10)))
  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
  }

  rbind(
    cbind(tau*model_data$lambdainv,-tau*model_data$lambdainv %*% Ad),
    cbind(-tau*t(Ad)%*% model_data$lambdainv,Suinv + tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,Ad)))
  )
}

Q_matrix_both_rw2_linear <- function(theta,model_data,tau = exp(7), buffer = 1/exp(7)) {

  if (is.null(model_data$A)) stop("no rw2 components in model")

  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)

  howmanyrw2 <- length(whichrw2)
  # The thetas are in order.


  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% map("Ad") %>% purrr::reduce(cbind)

  # APPLY ZERO CONSTRAINTS
  # Directly apply constraint that U_t = 0 by deleting the t^th column of A,
  # and the t^th row and column of Suinv

  if (!(0 %in% model_data$vectorofcolumnstoremove)) {
    Ad <- Ad[ ,-model_data$vectorofcolumnstoremove]
    Suinv <- Suinv[-model_data$vectorofcolumnstoremove,-model_data$vectorofcolumnstoremove]
  }

  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(10)))
  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
  }

  # Linear
  if (model_data$p == 0) stop("no linear terms in model")
  Sbinv <- Matrix::Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))

  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Ad
  Q13 <- -tau*model_data$lambdainv %*% model_data$Xd
  Q22 <- Suinv + tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,Ad))
  Q23 <- tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,model_data$Xd))
  Q33 <- Sbinv + tau*Matrix::crossprod(Matrix::crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd)
  rbind(
    cbind(Q11,Q12,Q13),
    cbind(t(Q12),Q22,Q23),
    cbind(t(Q13),t(Q23),Q33)
  )
}

Q_matrix_both_rw2_iid <- function(theta,model_data,tau = exp(7), buffer = 1/exp(7)) {


  if (is.null(model_data$A)) stop("no rw2 components in model")

  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)

  howmanyrw2 <- length(whichrw2)
  # The thetas are in order.


  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% map("Ad") %>% purrr::reduce(cbind)

  # APPLY ZERO CONSTRAINTS
  # Directly apply constraint that U_t = 0 by deleting the t^th column of A,
  # and the t^th row and column of Suinv

  if (!(0 %in% model_data$vectorofcolumnstoremove)) {
    Ad <- Ad[ ,-model_data$vectorofcolumnstoremove]
    Suinv <- Suinv[-model_data$vectorofcolumnstoremove,-model_data$vectorofcolumnstoremove]
  }

  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(12)))

  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
  }

  # IID
  if (is.null(model_data$B)) stop("no iid components in model")
  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)

  howmanyiid <- length(whichiid)
  # The thetas are in order.

  Suinv2 <- purrr::map2(whichiid,1:howmanyiid,
                        ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)

  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Ad
  Q13 <- -tau*model_data$lambdainv %*% Bd
  Q22 <- Suinv + tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,Ad))
  Q23 <- tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,Bd))
  Q33 <- Suinv2 + tau*Matrix::crossprod(Matrix::crossprod(model_data$lambdainv,Bd),Bd)
  rbind(
    cbind(Q11,Q12,Q13),
    cbind(t(Q12),Q22,Q23),
    cbind(t(Q13),t(Q23),Q33)
  )
}


Q_matrix_both_iid_linear <- function(theta,model_data,tau = exp(7)) {


  if (is.null(model_data$B)) stop("no iid components in model")
  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)

  howmanyiid <- length(whichiid)
  # The thetas are in order.

  Suinv2 <- purrr::map2(whichiid,1:howmanyiid,
                        ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)



  # linear
  if (model_data$p == 0) stop("no linear terms in model")
  Sbinv <- Matrix::Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))

  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Bd
  Q13 <- -tau*model_data$lambdainv %*% model_data$Xd
  Q22 <- Suinv2 + tau * Matrix::crossprod(Bd,Matrix::crossprod(model_data$lambdainv,Bd))
  Q23 <- tau * Matrix::crossprod(Bd,Matrix::crossprod(model_data$lambdainv,model_data$Xd))
  Q33 <- Sbinv + tau*Matrix::crossprod(Matrix::crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd)
  rbind(
    cbind(Q11,Q12,Q13),
    cbind(t(Q12),Q22,Q23),
    cbind(t(Q13),t(Q23),Q33)
  )
}


# Q matrix for model containing all of linear, rw2 and iid terms
Q_matrix_all <- function(theta,model_data,tau = exp(7)) {

  if (is.null(model_data$A)) stop("no rw2 components in model")

  whichrw2 <- model_data$modelspec %>% dplyr::filter(model == "rw2") %>% pull(covariate)

  howmanyrw2 <- length(whichrw2)
  # The thetas are in order.


  Suinv <- purrr::map2(whichrw2,1:howmanyrw2,
                       ~Q_matrix_rw2_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()


  # The full RE design matrix is the cbinded ones from each sub model
  Ad <- model_data$A %>% map("Ad") %>% purrr::reduce(cbind)

  # APPLY ZERO CONSTRAINTS
  # Directly apply constraint that U_t = 0 by deleting the t^th column of A,
  # and the t^th row and column of Suinv

  if (!(0 %in% model_data$vectorofcolumnstoremove)) {
    Ad <- Ad[ ,-model_data$vectorofcolumnstoremove]
    Suinv <- Suinv[-model_data$vectorofcolumnstoremove,-model_data$vectorofcolumnstoremove]
  }

  Qe <- eigen(Suinv,symmetric = T)
  index <- which(Qe$values <= (1/exp(12)))

  for (i in 1:length(index)) {
    Suinv <- Suinv + buffer * Qe$vectors[,index[i]] %*% t(Qe$vectors[,index[i]])
  }

  # Linear
  if (model_data$p == 0) stop("no linear terms in model")
  Sbinv <- Matrix::Diagonal(ncol(model_data$Xd),exp(model_data$beta_logprec))

  #IID random effect
  if (is.null(model_data$B)) stop("no iid components in model")
  whichiid <- model_data$modelspec %>% dplyr::filter(model == "iid") %>% pull(covariate)

  howmanyiid <- length(whichiid)
  # The thetas are in order.

  Suinv2 <- purrr::map2(whichiid,1:howmanyiid,
                        ~Q_matrix_iid_one_component(theta,model_data,covariate = .x)) %>%
    Matrix::bdiag()

  # The full RE design matrix is the cbinded ones from each sub model
  Bd <- model_data$B %>% map("Bd") %>% purrr::reduce(cbind)


  # Construct the matrix
  Q11 <- tau*model_data$lambdainv
  Q12 <- -tau*model_data$lambdainv %*% Ad
  Q13 <- -tau*model_data$lambdainv %*% Bd
  Q14 <- -tau*model_data$lambdainv %*% model_data$Xd
  Q22 <- Suinv + tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,Ad))
  Q23 <- tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,Bd))
  Q24 <- tau * Matrix::crossprod(Ad,Matrix::crossprod(model_data$lambdainv,model_data$Xd))
  Q33 <- Suinv2 + tau * Matrix::crossprod(Bd,Matrix::crossprod(model_data$lambdainv,Bd))
  Q34 <- tau * Matrix::crossprod(Bd,Matrix::crossprod(model_data$lambdainv,model_data$Xd))
  Q44 <- Sbinv + tau*Matrix::crossprod(Matrix::crossprod(model_data$lambdainv,model_data$Xd),model_data$Xd)
  rbind(
    cbind(Q11,Q12,Q13,Q14),
    cbind(t(Q12),Q22,Q23,Q24),
    cbind(t(Q13),t(Q23),Q33,Q34),
    cbind(t(Q14),t(Q24),t(Q34),Q44)
  )
}



# Make a general function to compute the Q matrix for a model
# This can be written to just call one of the above... somehow. Maybe as an option inside model_data
Q_matrix <- function(theta,model_data,tau = exp(7),forcesymm = TRUE) {
  # theta: vector of hyperparameters. The structure of this will depend on the model,
  # as specified by model_data
  if (is.null(model_data$A)) {
    if (model_data$p == 0) {stop("both X and A are null...")}
    else if(is.null(model_data$B)){
      mat <- Q_matrix_linear(theta,model_data,tau)
    }
    else{
      mat <- Q_matrix_both_iid_linear(theta,model_data,tau)
    }
  }
  else {
    if (model_data$p > 0) {
      if(is.null(model_data$B)){
        mat <- Q_matrix_both_rw2_linear(theta,model_data,tau)
      }
      else{
        mat <- Q_matrix_all(theta,model_data,tau)
      }
    }
    else {
      if(is.null(model_data$B)){
        mat <- Q_matrix_rw2(theta,model_data,tau)
      }
      # warning("random walk-only models are rank deficient for case crossover. adding fudge factor. consider adding a linear term.")
      else{
        mat <- Q_matrix_both_rw2_iid(theta,model_data,tau)
      }
      # mat <- mat + (1/tau) * Diagonal(dim(mat)[1])
    }
  }
  if (forcesymm) return(forceSymmetric(mat))
  mat
}






### Priors and Posteriors ###

#' log-prior of W
#'
#' @description These functions calculate priors and posteriors of the latent Gaussian variables.
#'
#' @param W Value of W = (delta,gamma,beta) to calculate the prior/posterior at. The order is as
#' it appears in the appendix of the paper- (delta_1_1,...,delta_n_Jn,gamma_1,...,gamma_M,beta_1,...,beta_p).
#' The order of gamma and beta is the same as they are listed in model_data$model_elements.
#' @param model_data A ccmodeldata object returned by model_setup().
#' @param theta Value of the hyperparameter vector theta. The Q matrix depends on this. Can leave as NULL
#' if you're passing in the Q matrix
#' @param Q The Q-matrix as returned by Q_matrix(model). If not provided, will be calculated.
#'
logprior_W <- function(W,theta,model_data,Q = NULL) {
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  -as.numeric((1/2)*crossprod(W,crossprod(Q,W)))
}


#' log-prior of W when we use AGHQ
#'
#' @description These functions calculate priors and posteriors of the latent Gaussian variables.
#'
#' @param W Value of W = (delta,gamma,beta) to calculate the prior/posterior at. The order is as
#' it appears in the appendix of the paper- (delta_1_1,...,delta_n_Jn,gamma_1,...,gamma_M,beta_1,...,beta_p).
#' The order of gamma and beta is the same as they are listed in model_data$model_elements.
#' @param model_data A ccmodeldata object returned by model_setup().
#' @param theta Value of the hyperparameter vector theta. The Q matrix depends on this. Can leave as NULL
#' if you're passing in the Q matrix
#' @param Q The Q-matrix as returned by Q_matrix(model). If not provided, will be calculated.
#'
logprior_W_aghq <- function(W,theta,model_data,Q = NULL) {
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  -as.numeric((1/2)*crossprod(W,crossprod(Q,W))) + 0.5 * as.numeric(determinant(Q,logarithm = TRUE)$modulus)
}

#' Log-posterior of W, given theta and y
#'
#' @rdname logprior_W
#'
log_posterior_W <- function(W,theta,model_data,Q = NULL) {
  logprior_W(W,theta,model_data,Q) + log_likelihood(W,model_data)
}


#' Gradient of log-posterior of W, given theta and y
#'
#' @rdname logprior_W
#' @param model_data Model's data created by abcoxph_setup
#'
grad_log_posterior_W <- function(W,theta,model_data,Q = NULL) {
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  -as.numeric(Matrix::crossprod(Q,W)) + grad_log_likelihood(W,model_data)
}

#' Hessian of log-posterior of W, given theta and y
#'
#' @rdname logprior_W
#' @param model_data Model's data created by abcoxph_setup
#'
hessian_log_posterior_W <- function(W,theta = NULL,Q = NULL, model_data) {
  if (is.null(theta) & is.null(Q)) stop("One of Q or theta must be provided")
  if (is.null(Q)) Q <- Q_matrix(theta,model_data)
  A <- as.matrix(-(Q + hessian_log_likelihood(W,model_data)))
  A <- as(A,'dgCMatrix')
  A
}




### HYPERPARAMETERS ###

#' Log-posterior approximation for theta and sigma
#'
#' @description Compute the log-posterior for theta, the log-precision, and sigma, the standard deviation, of the RW2 model components.
#'
#' @param theta Vector of log-precisions at which to evaluate the log-posterior. theta = -2*log(sigma).
#' @param W The mode of the log-posterior of W|theta,y for this theta. Computed outside and passed in.
#' @param model_data a ccmodeldata object created using model_setup()
#' @param Q Optional, Q matrix evaluated at this theta. Will be created if not supplied.
#'
log_posterior_theta <- function(theta,W,model_data,Q = NULL) {
  # W is the mode of log_posterior_W(theta)
  if (is.null(Q)) {
    Q <- Q_matrix(theta,model_data)
  }
  Q_p_C <- -1 * hessian_log_posterior_W(W,Q = Q,model_data = model_data)
  term1 <- log_likelihood(W,model_data)
  dt <- determinant(Q,logarithm = TRUE) # For this, we DO need the determinant (original)
  # term2_det <- (1/2) * as.numeric(dt$modulus * dt$sign)
  term2_det <- (1/2) * as.numeric(dt$modulus) #(original)
  term2 <- logprior_W(W,theta,model_data) # Doesn't contain the determinant
  term3 <- model_data$theta_logprior(theta)
  qcdet <- determinant(Q_p_C,logarithm = TRUE)
  # term4 <- -(1/2)*as.numeric(qcdet$modulus * qcdet$sign) # The gaussian approx evaluated at conditional mode
  term4 <- -(1/2) * as.numeric(qcdet$modulus) # (original)
  as.numeric(term1 + term2_det + term2 + term3 + term4)
}

#' Log-posterior approximation for theta and sigma
#'
#' @rdname log_posterior_theta
#' @param sigma Vector of standard deviations at which to evaluate the log-posterior. sigma = exp(-.5 * theta)
#'
log_posterior_sigma <- function(sigma,W,model_data,Q = NULL) {
  length(sigma)* log(2) - sum(log(sigma)) +
    log_posterior_theta(theta = -2 * log(sigma),
                        W = W,
                        model_data = model_data,
                        Q = Q)
}







