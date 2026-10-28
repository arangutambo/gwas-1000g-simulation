# REML estimate of SNP heritability from a GRM, using its eigendecomposition.
# Model: y = X b + g + e, g ~ N(0, s_g * A), e ~ N(0, s_e * I).
# In the eigenbasis of A the covariance is diagonal, so each likelihood evaluation is O(n p^2).
reml_h2 <- function(y, eig, X = NULL) {
  n <- length(y)
  X <- if (is.null(X)) matrix(1, n, 1) else cbind(1, as.matrix(X))
  ys <- crossprod(eig$U, y)
  Xs <- crossprod(eig$U, X)
  lam <- pmax(eig$lam, 0)
  nll <- function(par) {
    v <- exp(par[1]) * lam + exp(par[2])
    w <- 1 / v
    XtWX <- crossprod(Xs, Xs * w)
    beta <- solve(XtWX, crossprod(Xs, ys * w))
    r <- ys - Xs %*% beta
    0.5 * (sum(log(v)) + as.numeric(determinant(XtWX, logarithm = TRUE)$modulus) + sum(r^2 * w))
  }
  v0 <- var(y)
  o <- optim(log(c(v0 / 2, v0 / 2)), nll, method = "Nelder-Mead",
             control = list(reltol = 1e-10, maxit = 1000))
  sg <- exp(o$par[1]); se <- exp(o$par[2])
  c(sigma_g = sg, sigma_e = se, h2 = sg / (sg + se), converged = as.numeric(o$convergence == 0))
}
