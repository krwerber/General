#Writing functions for both fixed and random effects estimators

require(plm)
require(magrittr)
require(matlab)
require(dplyr)

within_estimator <- function(X, y) {
  if (class(X)[1] != "pdata.frame") {
    stop("X must be pdata.frame")
  }
  
  dims <- pdim(X)
  n <- dims$nT$n
  t <- dims$nT$T
  
  I_n <- diag(n)
  I_t <- diag(t)
  I_nt <- diag(n * t)
  l_n <- t(ones(1, n))
  l_t <- t(ones(1, t))
  l_nt <- t(ones(1, n * t))
  
  Z_mu <- I_n %x% l_t
  Z_lamda <- l_n %x% I_t
  
  P_mu <- I_n %x% ((1 / t) * l_t %*% t(l_t))
  P_lamda <- ((1 / n) * l_n %*% t(l_n)) %x% I_t
  P_mu_lamda <- ((1 / n) * l_n %*% t(l_n)) %x% ((1 / t) * l_t %*% t(l_t))
  #P_mu <- Z_mu %*% solve(t(Z_mu) %*% Z_mu) %*% t(Z_mu)
  #P_lamda <- Z_lamda %*% solve(t(Z_lamda) %*% Z_lamda) %*% t(Z_lamda)
  #P_mu_lamda <- P_mu %*% P_lamda
  
  Q <- I_nt - P_mu - P_lamda + P_mu_lamda
  
  beta_hat <- solve(t(X) %*% Q %*% as.matrix(X)) %*% t(X) %*% Q %*% as.matrix(y)
  colnames(beta_hat) <- "estimates"
  
  ids <- dims$panel.names$id.names
  times <- dims$panel.names$time.names
  
  mu_i <-  array(NA, dim = n)
  for (i in 1:n) {
    mu_i[i] = mean(y[STATE == ids[i], 1]) - (t(colMeans(X[STATE == ids[i], ])) %*% beta_hat)
  }
  
  lamda_t <- array(NA, dim = t)
  for (i in 1:t) {
    lamda_t[i] = mean(y[YR == times[i], 1]) - (t(colMeans(X[YR == times[i], ])) %*% beta_hat)
  }
  
  nu = matrix(NA, nrow = n, ncol = t)
  
  for (i in 1:n) {
    for (j in 1:t) {
      nu[i, j] <- y[STATE == ids[i] & YR == times[j], 1]  - lamda_t[j] - 
        (as.matrix(X[STATE == ids[i] & YR == times[j], ]) %*% beta_hat)[1] - mu_i[i]
    }
  }
  
  sigma_squared_hat <- t(nu) %*% nu / ((n - 1) * (t - 1) - ncol(X_select))
  
  return(list(coefs = beta_hat,
              variance = sigma_squared_hat))
}

munnell <- readr::read_csv("~/Desktop/Spring2020/econ387/munnell.csv", col_names = TRUE) %>%
  pdata.frame(x = ., index = c("STATE", "YR"))

attach(munnell)

X <- data.frame(STATE, YR, log(P_CAP + HWY + WATER + UTIL), log(PC), log(EMP), UNEMP) %>%
  pdata.frame(x = ., index = c("STATE", "YR"), drop.index = TRUE)

y <- data.frame(STATE, YR, log(GSP)) %>%
  pdata.frame(x = ., index = c("STATE", "YR"), drop.index = TRUE)

b <- within_estimator(X, y)

X <- data.frame(STATE, YR, log(P_CAP + HWY + WATER + UTIL), log(PC), log(EMP), UNEMP) %>%
  pdata.frame(x = ., index = c("STATE", "YR"))

yr_dummies <- model.matrix(~ -1 + X$YR)[, -1]
state_dummies <- model.matrix(~ -1 + X$STATE)[, -1]

X_dummy <- cbind(X, yr_dummies, state_dummies) %>%
  pdata.frame(x = ., index = c("STATE", "YR"), drop.index = TRUE)

c <- within_estimator(X_dummy, y)
#matrix was not invertible

oneway_ire <- function(X, y) {
  if (class(X)[1] != "pdata.frame") {
    stop("X must be pdata.frame")
  }
  
  dims <- pdim(X)
  n <- dims$nT$n
  t <- dims$nT$T
  
  I_n <- diag(n)
  I_t <- diag(t)
  l_t <- ones(t, 1)
  l_nt <- ones(n * t, 1)
  
  X_emph <- cbind(l_nt, as.matrix(select(X, -STATE, -YR)))
  
  J_t <- l_t %*% t(l_t)
  J_t_bar <- (1 / t) * J_t
  
  Z_mu <- I_n %x% l_t
  
  Q_t <- I_t - J_t_bar
  
  P_mu <- I_t %x% J_t_bar
  Q_mu <- I_n %x% Q_t
  
  log_likelihood <- function(theta, sigmasq_mu, sigmasq_nu) {
    (-n * t * 0.5) * log(2 * pi) - 0.5 * log(sigmasq_mu ^ (n * (t - 1)) * (t * sigmasq_mu + sigmasq_nu)^n) -
      0.5 * t(y - X_emph %*% theta) %*% ((t * sigmasq_mu + sigmasq_nu)^(-1) * P_mu + sigmasq_nu^(-1) * Q_mu) %*%
      (y - X_emph %*% theta)
  }
  
  optim(par = matrix(data = 1, nrow = n * t, ncol = nrow(P_mu)), sigmasq_mu = 1, sigmasq_nu = 1, fn = log_likelihood, method = "BFGS")
}

a2 <- oneway_ire(X, y)
