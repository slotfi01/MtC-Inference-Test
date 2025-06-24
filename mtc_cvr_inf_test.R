mtc_cvr_inf_test <- function(b, B, rt1, rt0, alpha) {
  library(MASS) # For mvrnorm
  library(lpSolve) # For linear programming
  
  # Inputs
  ret_matrx <- cbind(rt1, rt0)
  S <- nrow(ret_matrx)
  
  # Calculate CVaR and MTC
  linprog_solver <- function(rt, S, alpha) {
    f0 <- c(rep(1/S, S), (1 - alpha))
    A0 <- cbind(-diag(S), -1)
    b0 <- rt
    const.dir <- rep("<=", length(b0)) 
    A_lower <- diag(S + 1)  # Identity matrix to enforce x >= lower bounds
    b_lower <- c(rep(0, S), -Inf)  # Lower bounds for the variables
    dir_lower <- rep(">=", S + 1)  # Directions for lower bounds
    
    # Combine all constraints
    A_full <- rbind(A0, A_lower)
    b_full <- c(b0, b_lower)
    dir_full <- c(const.dir, dir_lower)
    res <- lp(
      direction = "min",
      objective.in = f0,
      const.mat = A0,
      const.dir = const.dir,
      const.rhs = b0
    )
    return(res$objval / (1 - alpha))
  }
  
  cvr_1 <- linprog_solver(rt1, S, alpha)
  cvr_0 <- linprog_solver(rt0, S, alpha)
  mtc_1 <- mean(rt1) / cvr_1
  mtc_0 <- mean(rt0) / cvr_0
  
  # Bootstrap sampling
  rsb0 <- S - b + 1
  rsp_all <- replicate(B, sample(1:rsb0, rsb0, replace = TRUE))
  k <- ceiling(S / b)
  rsp_k <- apply(rsp_all, 2, function(x) sample(x, k, replace = TRUE))
  
  rt1star <- matrix(NA, nrow = S, ncol = B)
  rt0star <- matrix(NA, nrow = S, ncol = B)
  for (m in 1:B) {
    rt1boot <- c()
    rt0boot <- c()
    for (i in 1:nrow(rsp_k)) {
      minv <- rsp_k[i, m]
      maxv <- rsp_k[i, m] + b - 1
      sample <- seq(minv, maxv)
      rt1boot <- c(rt1boot, ret_matrx[sample, 1])
      rt0boot <- c(rt0boot, ret_matrx[sample, 2])
    }
    rt1star[, m] <- rt1boot[1:S]
    rt0star[, m] <- rt0boot[1:S]
  }
  
  # Step 3: Calculate YbStar
  YbStar <- matrix(NA, nrow = B, ncol = 4)
  for (m in 1:B) {
    ret_1 <- rt1star[, m]
    ret_0 <- rt0star[, m]
    CVaR_0 <- linprog_solver(ret_0, S, alpha)
    CVaR_1 <- linprog_solver(ret_1, S, alpha)
    YbStar[m, ] <- c(mean(ret_1), CVaR_1, mean(ret_0), CVaR_0)
  }
  
  # Step 4: Calculate MTC statistics
  mtc_diff <- mtc_1 - mtc_0
  mtcdiff_s <- sqrt(S) * mtc_diff
  
  Mbar <- colMeans(YbStar)
  SigmaStar_r <- cov(YbStar)
  
  c_und <- c(1/cvr_1, -0.5 * ((mtc_1 + mtc_0) / cvr_1), -1/cvr_0, 0.5 * ((mtc_1 + mtc_0) / cvr_0))
  tau2 <- S * t(c_und) %*% SigmaStar_r %*% c_und
  
  # Calculate p-value
  if (abs(mtc_diff) <= 1e-6) {
    mtc_pval <- 1
  } else {
    mtc_pval <- 1 - pnorm(mtcdiff_s, mean = 0, sd = sqrt(tau2))
  }
  
  # Step 4 for CVaR
  cvr_diff <- cvr_1 - cvr_0
  cvrdiff_s <- sqrt(S) * cvr_diff
  
  e_und <- c(0, 1, 0, -1)
  nu2 <- S * t(e_und) %*% SigmaStar_r %*% e_und
  
  # Calculate CVaR p-value
  if (abs(cvr_diff) <= 1e-6) {
    cvr_pval <- 1
  } else {
    cvr_pval <- 1 - pnorm(cvrdiff_s, mean = 0, sd = sqrt(nu2))
  }
  
  # Return results
  return(list(mtc_diff = mtc_diff, mtc_pval = mtc_pval, cvr_diff = cvr_diff, cvr_pval = cvr_pval))
}
