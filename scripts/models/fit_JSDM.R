#+ setup, message=FALSE
library(arm)
library(jagsUI)
library(Rcpp)

#+ data
data <- list(
  Y = y, # Matrix: columns = taxa, rows = sites
  X = cbind(1, X), # Matrix: column 1 = 1s, column 2--N = scaled site variables e.g., cbind(1, scale(X))
  covx = cov(X), # Covariance matrix of X, e.g., cov(X)
  K = 1 + n_env_vars, # Integer: 1 + number of site variables
  J = n_species, # Integer: number of taxa
  n = nrow(X), # Integer: number of sites
  I = diag(n_species), # Identity matrix, e.g. diag(J),
  df = n_species + 1 # Integer, Degrees of freedom e.g., J + 1
)

#+ model
run_model <- function(data, n.chains = 1, n.iter = 1000, n.adapt = 500, n.thin = 1) {
  jsdm_jags <- function() {
    model.file <- tempfile()
    cat(
      "model {
      for (i in 1:n) {
      Z[i, 1:J] ~ dmnorm(Mu[i, ], Tau)
      for (j in 1:J) {
      Mu[i, j] <- inprod(B_raw[j, ], X[i, ])
      Y[i, j] ~ dbern(step(Z[i, j]))
      }
      }
      for (j in 1:J) {
      sigma_[j] <- sqrt(Sigma[j, j])
      env_sigma_[j] <- sqrt(EnvSigma[j, j])
      for (k in 1:K) {
      B_raw[j, k] ~ dnorm(mu[k], tau[k])
      B[j, k] <- B_raw[j, k] / sigma_[j]
      }
      for (j_ in 1:J) {
      Rho[j, j_] <- Sigma[j, j_] / (sigma_[j] * sigma_[j_])
      EnvRho[j, j_] <- EnvSigma[j, j_] / (env_sigma_[j] * env_sigma_[j_])
      EnvSigma[j, j_] <-
      sum(EnvSigma1[, j, j_]) + sum(EnvSigma2[, , j, j_])
      for (k in 2:K) {
      EnvSigma1[k - 1, j, j_] <- B[j, k] * B[j_, k]
      for (k_ in 2:K) {
      EnvSigma2[k - 1, k_ - 1, j, j_] <-
      B[j, k] * B[j_, k_] * ifelse(k_ != k, covx[k, k_], 0)
      }
      }
      }
      }
      for (k in 1:K) {
      mu[k] ~ dnorm(0, 1)
      tau[k] <- pow(sigma[k], -2)
      sigma[k] ~ dnorm(0, 1)T(0,)
      }
      Tau ~ dwish(I, df)
      Sigma <- inverse(Tau)
  }",
      file = model.file
    )
    model.file
    }
  
  inits <- function(data) {
    Y <- as.matrix(data$Y)
    X <- as.matrix(data$X)[, -1]
    Tau <- rWishart(1, data$df, data$I)[, , 1]
    Sigma <- solve(Tau)
    Z <- rep(0, data$J)
    Z <- mvrnorm(1, Z, Sigma)
    Z <- replicate(data$n, Z)
    Z <- t(Z)
    Z <- abs(Z)
    Z <- ifelse(Y, Z, -Z)
    Sigma <- cov(Z)
    B <- sapply(
      seq_len(data$J),
      function(x) coef(bayesglm(Y[, x] ~ X, family = binomial(link = "probit")))
    )
    B <- t(B)
    B_raw <- B * sqrt(diag(Sigma))
    mu <- apply(B_raw, 2, mean)
    sigma <- pmin(99, apply(B_raw, 2, sd))
    Tau <- solve(Sigma)
    list(Tau = Tau, Z = Z, B_raw = B_raw, mu = mu, sigma = sigma)
  }
  
  jags(
    data,
    function() inits(data), c("B", "Rho", "EnvRho"), jsdm_jags(),
    n.chains = n.chains, n.iter = n.iter, n.adapt = n.adapt, n.thin = n.thin,
    parallel = TRUE, DIC = FALSE
  )
  }

JSDM <- run_model(data,
                  n.chains = n.chains,
                  n.iter = n.iter,
                  n.adapt = n.burn,
                  n.thin = n.thin)
