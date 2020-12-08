# This file contains miscillaneous top-level non-exported functions
# for the abcoxph package


# INTERNAL: parse formula
parse_formula <- function(ff) {
  # Parse the formula ff into linear and smooth terms and a id
  # Linear terms will be passed to model.matrix()
  # Smooth terms will be handled in a proprietary manner
  fullff <- ff # Save the input formula for later printing.

  # Grab the RHS elements of the formula
  ff_elements <- attributes(terms(ff))$term.labels
  # Get the names of the variables. This strips off the s() wrappers
  # and includes the respose.
  ff_variables <- all.vars(ff)

  # Split the names into response, linear, and smooth terms.
  # Throw an error if any terms are not supported.
  response <- ff_variables[1] # Response is always the first
  # Match smooth terms using regex
  # Will only match terms that look like s(variable). Extra stuff in the s()
  # will cause an error. Keep formulas clean!
  smooth <- stringr::str_extract(ff_elements,"^s\\(\\w+\\)$")
  smooth <- smooth[!is.na(smooth)]
  # Remove it from ff
  if (length(smooth) > 0) {
    for (tm in smooth) ff <- update(ff,formula(stringr::str_c(".~.-",tm)))
  }
  # Strip off the s() part
  smooth <- stringr::str_remove(smooth,"^s\\(") %>% stringr::str_remove("\\)$")
  # Match id terms in same way as smooth
  id <- stringr::str_extract(ff_elements,"^id\\(\\w+\\)$")
  id <- id[!is.na(id)]
  # Remove it from ff
  if (length(id) > 0) {
    for (tm in id) ff <- update(ff,formula(stringr::str_c(".~.-",tm)))
  }
  # Strip off the s() part
  id <- stringr::str_remove(id,"^id\\(") %>% stringr::str_remove("\\)$")

  # All the terms that are left are treated as linear and passed to model.matrix
  # within the model_setup function. Any errors will come from that call.

  # Finally, remove the intercept. If a -1 was already in the formula, this operation does nothing,
  # so it's safe.
  ff <- update(ff,.~.-1)

  # Return a named list of vectors with the names of each type of term
  list(
    # linear = attributes(terms(ff))$term.labels,
    linear = setdiff(all.vars(ff),response),
    linear_formula = ff,
    smooth = smooth,
    id = id,
    response = response,
    call = fullff
  )
}

# INTERNAL: get the degree of a polynomial from a formula
get_polynomial_degree <- function(ff) {
  ffvars <- all.vars(ff)[-1]
  ffattr <- attributes(terms(ff))$term.labels

  varnameregex <- "A-Za-z0-9_."

  degree_1 <- stringr::str_extract(ffattr,stringr::str_c("^[",varnameregex,"]+$"))
  degree_1 <- degree_1[!is.na(degree_1)]

  degree_more_than_1 <- stringr::str_extract(ffattr,stringr::str_c("^poly\\([",varnameregex,"=\\s\\,\\)]*"))
  degree_more_than_1 <- degree_more_than_1[!is.na(degree_more_than_1)]

  # Get the names
  deg_mt1_names <- stringr::str_extract(degree_more_than_1,stringr::str_c("^poly\\([",varnameregex,"]+")) %>%
    stringr::str_remove("^poly\\(")

  deg_mt1_degrees <- stringr::str_extract(degree_more_than_1,stringr::str_c("^poly\\([",varnameregex,"]+\\,\\s?[A-Za-z\\s=]*[0-9]")) %>%
    stringr::str_remove(stringr::str_c("^poly\\([",varnameregex,"]+\\,\\s?[A-Za-z\\s=]*")) %>%
      as.numeric()

  out <- c(rep(1,length(degree_1)),deg_mt1_degrees)
  names(out) <- c(degree_1,deg_mt1_names)
  out
}


# INTERNAL: prescribe default control arguments
abcox_default_control <- function() {
  list(
    smooth_prior = list(),
    linear_constraints = list(),
    doparallel = TRUE,
    thetaaccuracy = 3,
    sparsetheta = FALSE,
    thetarange = c(-1,1),
    beta_prior_logprec = log(1/10),
    opt_control = list(
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
  )
}



# INTERNAL: take a vector of values, a vector of indices, and return a vector
# containing the values of the first with zeroes spliced in at the indices
# specified by the second
stitch_zero_vector <- function(x,z) {
  if (any(z > length(x) + 1)) stop("You are asking for zeroes at positions not covered by x")
  y <- as(x,"sparseVector")
  # This drops zero values in x automatically- add them back in
  y@i <- 1:length(x)
  y@x <- x
  for (j in z) {
    y@i[y@i >= j] <- y@i[y@i >= j] + 1
  }
  as.numeric(y)
}

# INTERNAL: compute quantiles given a vector of grid points and logged density evaluations
compute_quantiles <- function(tt,pp,origgrid) {
  pp <- exp(pp)
}


### MODEL SETUP ###
# Function to create the difference matrix D
create_diff_matrix <- function(n) {
  # n is the total # sample size.
  cbind(Matrix::Matrix(1,n-1,1),Matrix::Diagonal(n-1,-1))
}
# Function to create the crossproduct's inverse,(DD^T)^(-1)
create_full_dtcp_matrix <- function(n) {
  m <- Matrix::Diagonal(n-1,1) - Matrix::Matrix(1/n,n-1,n-1)
  m
}

# Use to order your data first!!! It orders your dataset based on observed times:
arrange_data <- function(data){
  newdata <- dplyr::arrange(data,times)
  newdata
}

#Function that creates the adjust rank that will be used later for Breslow's adjustment:
Get_Adj <- function(model_data){
  new_data <- model_data
  new_data$adjRank <- rank(model_data$times,ties.method = "min")
  new_data
}

# Function to take a covariate and return the appropriate element of "Alist"(rw2)
create_alist_element <- function(u,constraint = NULL) {
  # u: covariate. NOT sorted and MAY contain ties/repeated values, in general.
  # constraint: vector containing values of u for which random effect U should be
  # constrained to be zero.
  lu <- length(u)
  A <- Matrix::Diagonal(n = lu)[match(u,unique(u)),order(unique(u))]
  model <- "rw2"
  constrzero <- NULL
  if (!is.null(constraint)) {
    constrzero <- match(constraint,sort(unique(u)))
    if (any(is.na(constrzero))) warning(paste0("no match found for constraint: ",constraint[which(is.na(constrzero))]))
  }
  list(u = u,A = A,model = model,constrzero = constrzero)
}

# Function to take a covariate and return the appropriate element of "Blist"(iid)
create_blist_element <- function(u,constraint = NULL) {
  # u: covariate. NOT sorted and MAY contain ties/repeated values, in general.
  # constraint: vector containing values of u for which random effect U should be
  # constrained to be zero.
  lu <- length(u)
  B <- Matrix::Diagonal(n = lu)[match(u,unique(u)),order(unique(u))]
  model <- "iid"
  constrzero <- NULL
  if (!is.null(constraint)) {
    constrzero <- match(constraint,sort(unique(u)))
    if (any(is.na(constrzero))) warning(paste0("no match found for constraint: ",constraint[which(is.na(constrzero))]))
  }
  list(u = u,B = B,model = model,constrzero = constrzero)
}




# Bin the covariates into 100 bins
# This gets the dimension of the RW part of the latent field down
# from 3919 + 3919 + 3922 = 11,760 to 100 x 3 = 300 (!)
bin_covariate <- function(u,bins = 100,type = "quantile",custombins = NULL,decimals = 5) {
  if (min(u) < 0) {
    lowerend <- min(u)*1.01
  } else {
    lowerend <- min(u) * .99
  }
  if (type == "quantile") {
    # bin into quantiles. quantile gives the upper end of the range.
    bininfo <- tibble::tibble(
      upper = stats::quantile(u,probs = (1:bins)/bins),
      lower = lag(upper,n = 1L,default = lowerend),
      midpoint = (lower + upper)/2
    )
  } else if (type == "equal") {
    # Bin into equally spaced bins
    bininfo <- tibble::tibble(
      upper = seq(min(u),max(u),length.out = bins),
      lower = lag(upper,n = 1L,default = lowerend),
      midpoint = (lower + upper)/2
    )
  } else if (type == "custom") {
    if (is.null(custombins)) stop("type == custom but no custom bins provided")
    # The provided bins are treated as upper end points
    # Make sure to provide a big range that covers all the data.
    lowerend <- min(custombins)
    bininfo <- tibble::tibble(
      upper = custombins,
      lower = lag(upper,n = 1L,default = lowerend),
      midpoint = (lower + upper)/2
    )
  }

  # Now if u[i] is between lower and upper, out[i] = midpoint
  out <- numeric(length(u))
  for (i in 1:length(u)) {
    tmp <- bininfo %>%
      filter(lower < u[i],upper >= u[i]) %>%
      dplyr::pull(midpoint)
    # if (length(tmp) == 0) cat(i,"\n")
    out[i] <- tmp
  }
  round(out,decimals)
}




