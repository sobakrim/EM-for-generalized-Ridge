library(Matrix)
library(geoR)
library(MASS)
library(corpcor)
library(mvnfast)

adjacency_matrix <- function(D, nv = 1) {
  # Function to create an adjacency matrix based on the distance matrix D
  # Args:
  #   D: A distance matrix (square matrix)
  #   nv: Number of nearest neighbors to connect (default is 1)
  # Returns:
  #   A sparse adjacency matrix
  
  # Initialize an empty adjacency matrix with the same dimensions as D
  adj_matrix <- matrix(0, nrow = nrow(D), ncol = ncol(D))
  
  # Loop through each row in the distance matrix
  for (i in 1:nrow(D)) {
    # Get the indices that would sort the row, excluding the first one (which is the point itself)
    nearest_indices <- order(D[i,])[-1]
    
    # Connect the current point to its nv nearest neighbors
    adj_matrix[i, nearest_indices[1:nv]] <- 1
    adj_matrix[nearest_indices[1:nv], i] <- 1
  }
  
  # Return the adjacency matrix as a sparse matrix
  return(Matrix::Matrix(adj_matrix, sparse = TRUE))
}

matern_cov <- function(par, D) {
  # Function to compute the Matern covariance matrix
  # Args:
  #   par: A vector of parameters for the Matern function (phi and kappa)
  #   D: A distance matrix (square matrix)
  # Returns:
  #   A Matern covariance matrix
  
  return(geoR::matern(D, phi = par[1], kappa = abs(par[2])))
}
precisCAR <- function(par, D) {
  # Function to compute the precision matrix for the Conditional Autoregressive (CAR) model
  # Args:
  #   par: A parameter for the CAR model
  #   D: A distance matrix (square matrix)
  # Returns:
  #   The precision matrix for the CAR model
  
  w <- adjacency_matrix(D, nv = 4)
  diag_inv_w <- diag(1 / rowSums(w))
  return(Matrix(diag(rowSums(w)) %*% (diag(1, ncol(D)) - par * diag_inv_w %*% w), sparse = TRUE))
}

simulate_beta <- function(D, par, mu = NULL, x = NULL, cov.model = c("matern", "CAR", "ridge")) {
  # Function to simulate beta values based on a covariance model
  # Args:
  #   D: A distance matrix (square matrix)
  #   par: A vector of parameters for the covariance function
  #   mu: Mean vector (default is NULL)
  #   x: Not used in this function (default is NULL)
  #   cov.model: Covariance model to use ("matern", "CAR", or "ridge")
  # Returns:
  #   A simulated vector of beta values
  
  require(geoR)
  
  # Determine the covariance matrix based on the selected model
  if (cov.model == "ridge") {
    sigma <- diag(par[1], ncol(D))
  } else if (cov.model == "matern") {
    sigma <- par[1] * matern_cov(par[2:3], D)
  } else if (cov.model == "CAR") {
    sigma <- par[1] * as.matrix(solve(precisCAR(par[2], D)))
  }
  
  # Simulate and return the beta values using the specified mean and covariance matrix
  return(MASS::mvrnorm(1, mu = mu, Sigma = sigma))
}

pseudo_det <- function(x) {
  # Function to compute the pseudo-determinant of a matrix
  # Args:
  #   x: A matrix
  # Returns:
  #   The pseudo-determinant of the matrix
  
  sv <- try(svd(x), silent = TRUE)
  if (!is.character(sv)) {
    sigma <- sv$d
    return(prod(sigma[sigma != 0]))
  } else {
    return(0)
  }
}

pseudo_inv <- function(sv) {
  # Function to compute the pseudo-inverse of a matrix from its SVD
  # Args:
  #   sv: A singular value decomposition list (containing u, d, and v)
  # Returns:
  #   The pseudo-inverse of the matrix
  
  d <- sv$d
  d[d != 0] <- 1 / d[d != 0]
  return(sv$v %*% diag(d) %*% t(sv$u))
}

loglik_matern <- function(par, x, D) {
  # Function to compute the log-likelihood for the Matern covariance model
  # Args:
  #   par: A vector of parameters for the Matern function
  #   x: Data vector
  #   D: Distance matrix (square matrix)
  # Returns:
  #   The log-likelihood value
  
  cor <- matern_cov(c(par, 3/2), D)
  if (sum(is.infinite(cor)) == 0 & sum(is.na(cor)) == 0) {
    sv <- svd(cor)
    inv <- pseudo_inv(sv)
    de <- as.numeric(determinant(inv)$modulus)
    s <- t(x) %*% inv %*% x
    if (!de == Inf & !de == -Inf & !is.na(de) & s > 0) {
      return(ncol(D) * log(s) - de)
    } else {
      return(1e+10 * abs(rnorm(1)))
    }
  } else {
    return(1e+10 * abs(rnorm(1)))
  }
}

loglik_CAR <- function(par, x, D) {
  # Function to compute the log-likelihood for the CAR model
  # Args:
  #   par: A vector of parameters for the CAR function
  #   x: Data vector
  #   D: Distance matrix (square matrix)
  # Returns:
  #   The log-likelihood value
  
  I <- precisCAR(par, D)
  if (sum(is.infinite(I)) == 0 & sum(is.na(I)) == 0) {
    de <- as.numeric(determinant(I)$modulus)
    s <- as.matrix(t(x) %*% I %*% x)
    if (!de == Inf & !de == -Inf & !is.na(de) & s > 0) {
      return(ncol(D) * log(s) - de)
    } else {
      return(1e+10 * abs(rnorm(1)))
    }
  } else {
    return(1e+10 * abs(rnorm(1)))
  }
}

simulate_all <- function(n, lonlat, par_beta, cov.model = "matern", par_x, sigma) {
  # Function to simulate data based on spatial covariates
  # Args:
  #   n: Number of observations
  #   lonlat: Matrix of longitude and latitude coordinates
  #   par_beta: Parameters for beta simulation
  #   cov.model: Covariance model to use ("matern" or "CAR")
  #   par_x: Parameters for x simulation
  #   sigma: Standard deviation for noise
  # Returns:
  #   A list containing simulated x, y, and beta
  
  D <- as.matrix(dist(lonlat))
  sigma_x <- par_x[1] * matern_cov(par_x[2:3], D)
  x <- rmvn(n, mu = rep(0, nrow(lonlat)), sigma = sigma_x)
  beta <- simulate_beta(D, par_beta, mu = rep(0, nrow(lonlat)), cov.model = cov.model)
  y <- x %*% beta + rnorm(n, sd = sigma)
  return(list(x = x, y = y, beta = beta))
}

init_em <- function(x, D, init, method = "gauss") {
  # Function to initialize the EM algorithm for parameter estimation
  # Args:
  #   x: Data vector
  #   D: Distance matrix (square matrix)
  #   init: Initial parameters for optimization
  #   method: Method for parameter estimation ("matern" or "CAR")
  # Returns:
  #   Estimated parameters
  
  if (method == "matern") {
    par <- c()
    par[2] <- optim(init, f = loglik_matern, x = x, D = D, control = list(maxit = 50))$par
    par[1] <- c((1/length(x)) * t(x) %*% corpcor::pseudoinverse(matern_cov(c(par[2], 5/2), D)) %*% x)
  }
  if (method == "CAR") {
    par <- c()
    par <- optim(init, f = loglik_CAR, x = x, D = D, method = "L-BFGS-B", lower = c(0, 0), upper = c(10, 1), control = list(maxit = 60))$par
  }
  return(par)
}

LogLikEM <- function(par, ebeta, D = NULL, cov.model = c("matern", "CAR")) {
  # Function to compute the log-likelihood for the EM algorithm
  # Args:
  #   par: A vector of parameters for the covariance function
  #   ebeta: Expected beta values
  #   D: Distance matrix (square matrix)
  #   cov.model: Covariance model to use ("matern" or "CAR")
  # Returns:
  #   The log-likelihood value
  
  if (cov.model == "matern") {
    cor <- matern_cov(c(par, 3/2), D)
    if (sum(is.infinite(cor)) == 0 & sum(is.na(cor)) == 0) {
      I <- try(chol(cor), silent = TRUE)
      if (!is.character(I[1])) {
        I <- chol2inv(I)
        de <- as.numeric(determinant(I)$modulus)
        s <- sum(diag(I %*% ebeta))
        if (!de == Inf & !de == -Inf & !is.na(de) & s > 0) {
          return(ncol(D) * log(s) - de)
        } else {
          return(10^10 * abs(rnorm(1)))
        }
      } else {
        return(10^10 * abs(rnorm(1)))
      }
    } else {
      return(10^10 * abs(rnorm(1)))
    }
  }
  
  if (cov.model == "CAR") {
    I <- precisCAR(par, D)
    if (sum(is.infinite(I)) == 0 & sum(is.na(I)) == 0) {
      de <- as.numeric(determinant(I)$modulus)
      A <- try(chol(I), silent = TRUE)
      s <- sum(diag(I %*% ebeta))
      if (!de == -Inf & !de == Inf & !is.na(de) & !is.character(de) & !is.character(A[1]) & s > 0) {
        return(ncol(D) * log(s) - de)
      } else {
        return(10^10 * abs(rnorm(1)))
      }
    } else {
      return(10^10 * abs(rnorm(1)))
    }
  }
}


em_gr <- function(x, y, folds = NULL, init = NULL, x_test = NULL, y_test = NULL, lonlat = NULL, optim_maxit = 100, cov.model = c("matern", "CAR", "ridge"), return_Likelihood = TRUE, maxit = 10, tol = 1e-1) {
  # Function to fit the EM algorithm for Generalized Ridge regression with spatial covariates
  # Args:
  #   x: Predictor matrix
  #   y: Response vector
  #   folds: Cross-validation folds (default is NULL)
  #   init: Initial parameters (default is NULL)
  #   x_test: Test predictor matrix (default is NULL)
  #   y_test: Test response vector (default is NULL)
  #   lonlat: Matrix of longitude and latitude coordinates
  #   optim_maxit: Maximum iterations for optimization (default is 100)
  #   cov.model: Covariance model to use ("matern", "CAR", or "ridge")
  #   return_Likelihood: Boolean indicating whether to return the likelihood (default is TRUE)
  #   maxit: Maximum iterations for the EM algorithm (default is 10)
  #   tol: Tolerance for convergence (default is 1e-1)
  # Returns:
  #   A list containing sigma_beta_inv, sigma, mu, LL, and cov_par
  
  LL <- c()
  xtx <- crossprod(x)
  n <- nrow(x)
  D <- as.matrix(dist(lonlat))
  
  if (cov.model == "CAR") {
    par <- c(0.1, 0.3)
    sigma_beta_inv <- (1 / par[1]) * precisCAR(par[2], D)
    sigma <- 1
  }
  
  if (cov.model == "matern") {
    if (is.null(init)) {
      par <- c(1, 10)
      sigma_beta_inv <- (1 / par[1]) * corpcor::pseudoinverse(matern_cov(c(par[2], 3/2), D))
      sigma <- 1
      sigma_beta_y <- corpcor::pseudoinverse(sigma_beta_inv + xtx / sigma)
      mu <- (1 / sigma) * sigma_beta_y %*% t(x) %*% y
      Ebeta <- sigma_beta_y + mu %*% t(mu)
      par[1] <- c((1 / ncol(D)) * sum(diag(as.matrix(sigma_beta_inv %*% Ebeta))))
      sigma_beta_inv <- (1 / par[1]) * corpcor::pseudoinverse(matern_cov(c(par[2], 3/2), D))
      sigma <- c((1 / n) * (t(y) %*% y - 2 * t(y) %*% x %*% mu + sum(diag(xtx %*% Ebeta))))
    } else {
      par <- init$par
      sigma_beta_inv <- (1 / par[1]) * corpcor::pseudoinverse(matern_cov(par[2:3], D))
      sigma <- init$sigma
    }
  }
  
  if (cov.model == "ridge") {
    sigma_beta_inv <- diag(0, ncol(x))
    b <- corpcor::pseudoinverse(xtx + diag(0, ncol(x))) %*% t(x) %*% y
    sigma <- 1
  }
  
  for (i in 1:maxit) {
    sigma_beta_y <- corpcor::pseudoinverse(as.matrix(sigma_beta_inv) + xtx / sigma)
    mu <- (1 / sigma) * sigma_beta_y %*% t(x) %*% y
    Ebeta <- sigma_beta_y + mu %*% t(mu)
    sigma <- c((1 / n) * (t(y) %*% y - 2 * t(y) %*% x %*% mu + sum(xtx * Ebeta)))
    
    if (cov.model == "ridge") {
      par <- ncol(x) / sum(diag(Ebeta))
      sigma_beta_inv <- diag(par, ncol(x))
    } else {
      if (cov.model == "CAR") {
        print(paste0("--------iteration ", i, " ------------"))
        par[2] <- optim(par = par[2], fn = LogLikEM, D = D, ebeta = Ebeta, cov.model = cov.model, method = "L-BFGS-B", lower = 0, upper = 1, control = list(maxit = optim_maxit, trace = 2))$par
        sigma_beta_inv <- precisCAR(par[2], D)
        par[1] <- c((1 / ncol(D)) * sum(diag(sigma_beta_inv %*% Ebeta)))
        sigma_beta_inv <- (1 / par[1]) * sigma_beta_inv
        print(par)
        Sys.sleep(1)
      }
      if (cov.model == "matern") {
        print(paste0("--------iteration ", i, " ------------"))
        par[2] <- optim(par = par[2], fn = LogLikEM, D = D, ebeta = Ebeta, cov.model = cov.model, method = "L-BFGS-B", lower = 1e-5, upper = 100, control = list(maxit = optim_maxit, trace = 2))$par
        sigma_beta_inv <- corpcor::pseudoinverse(matern_cov(c(par[2], 3/2), D))
        par[1] <- c((1 / ncol(D)) * sum(sigma_beta_inv * Ebeta))
        sigma_beta_inv <- (1 / par[1]) * sigma_beta_inv
        print(par)
        Sys.sleep(1)
      }
    }
    
    if (return_Likelihood) {
      if (cov.model == "matern") {
        sigmaa <- par[1] * matern_cov(c(par[2], 3/2), D)
        sigma_y <- x %*% sigmaa %*% t(x) + diag(sigma, length(y))
      }
      if (cov.model == "ridge") sigma_y <- x %*% diag(1 / par, ncol(x)) %*% t(x) + diag(sigma, length(y))
      if (cov.model == "CAR") sigma_y <- x %*% solve(sigma_beta_inv) %*% t(x) + diag(sigma, length(y))
      
      LL <- c(LL, dmvn(t(y), mu = rep(0, length(y)), sigma = as.matrix(sigma_y), log = TRUE))
    }
    
    if (return_Likelihood & i > 2) {
      if ((abs(LL[i] - LL[i - 1]) / ((abs(LL[i]) + abs(LL[i - 1]) + tol) / 2)) < tol) break
    }
  }
  
  return(list(sigma_beta_inv = sigma_beta_inv, sigma = sigma, mu = mu, LL = LL, cov_par = par))
}
