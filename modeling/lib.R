# For space-filling designs, modeling is done in the full space.
library(nloptr)
Rcpp::sourceCpp("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/modeling/MIM_GP.cpp")

# MIM kernel
Fit.MIM.alpha <- function(X, y) {
  p <- ncol(X)
  
  alpha <- rep(1,p)
  theta <- rep(1, p)
  ini <- c(alpha, theta)
  lb <- c(rep(1e-1, p), rep(1e-2, p))
  ub <- c(rep(1e+1, p), rep(1e+2, p))
  
  param <- nloptr(ini, nllh_fun_alpha, lb = lb, ub = ub,
                  opts = list(algorithm = "NLOPT_LN_BOBYQA",
                              xtol_rel = 1e-4,
                              print_level = 0),
                  X = X, y = y, p = p)
  alpha = param$solution[1:p]
  theta = param$solution[(p+1):(2*p)]
  ans <- list(alpha = alpha, theta = theta, nllh = param$objective)
  return(ans)
}

Fit.MIM.beta <- function(X, y, delta, t) {
  p <- ncol(X)
  
  beta <- 0.5/max(t)
  theta <- rep(1, p)
  ini <- c(beta, theta)
  lb <- c(0.05, rep(1e-2, p))
  ub <- c(1/max(t)-1e-4, rep(1e+2, p))
  
  param <- nloptr(ini, nllh_fun, lb = lb, ub = ub,
                  opts = list(algorithm = "NLOPT_LN_BOBYQA",
                              xtol_rel = 1e-4,
                              print_level = 0),
                  X = X, y = y, p = p,
                  delta = delta, t = t)
  beta <- param$solution[1]
  theta = param$solution[2:(p+1)]
  alpha <- -log(1 - beta*t)/log(1 + delta^2/theta^2)
  ans <- list(alpha = alpha, theta = theta, beta = beta, nllh = param$objective)
  return(ans)
}

# surrogate model
surrogate <- function(X, X_type, l, fun, seed, kernel, test_X) {
  set.seed(seed)
  p <- ncol(X)
  rdx <- sample(p,p) # randomly arrange columns to factors
  X_rdx <- X[,rdx]
  y <- fun(X_rdx)
  test_y = fun(test_X)
  screening_level <- 1e-3

  if(X_type == "OFAT") { # For OFAT designs
    X_rdx_OFAT <- X_rdx[1:(l*(p+1)),]
    y_OFAT <- y[1:(l*(p+1))]
    t = MOFAT::measure(X_rdx_OFAT, y_OFAT)$t
  } else { # For other designs
    screen_fit = rkriging::Fit.Kriging(X_rdx, y, kernel.parameters=list(type="Gaussian"))
    m = 10000
    A=matrix(runif(m*p),nrow=m)
    B=matrix(runif(m*p),nrow=m)
    colnames(A) <- paste0("X", 1:p)
    colnames(B) <- paste0("X", 1:p)
    pred_fun <- function(new_X){rkriging::Predict.Kriging(screen_fit, new_X)$mean}
    a.sen=sensitivity::soboljansen(model = pred_fun, X1 = A, X2 = B)
    t=a.sen$T[,1]
  }

  norm_t = t / sum(t)
  active_col = which(norm_t > screening_level)
  active_norm_t = norm_t[active_col]
  new_p = length(active_col)
  active_X = matrix(X_rdx[,active_col], ncol = new_p)
  active_X_row = which(!duplicated(active_X)) # unique row idx in active_X
  active_X = active_X[active_X_row,,drop = FALSE]
  y = y[active_X_row]
  test_X = test_X[, active_col, drop = FALSE]

  if(kernel == "MIM") {
    if(X_type == "OFAT") {
      delta <- abs(X[1,] - X[l+1,])[which(abs(X[1,] - X[l+1,]) != 0)]
      fit <- Fit.MIM.beta(active_X, y, delta, active_norm_t)
    } else {
      fit <- Fit.MIM.alpha(active_X, y)
    }
    pred <- gp_predict(active_X, y, test_X, fit$alpha, fit$theta)
  } else if(kernel == "Gaussian") {
    fit <- rkriging::Fit.Kriging(active_X, y, kernel.parameters = list(type = "Gaussian"))
    pred <- rkriging::Predict.Kriging(fit, test_X)
  }
  mse <- mean((test_y - pred$mean)^2)

  res = list(mse = mse, nmse = mse / var(y), norm_t = norm_t, active_col = active_col)
  return(res)
}


### parallel for surrogate
library(future)
library(future.apply)
library(progressr)
workers <- max(1L, future::availableCores() - 2L)
future::plan(future::multisession, workers = workers)
# future::plan(future::multicore, workers = workers)
progressr::handlers("txtprogressbar")
cpp_file <- "/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/modeling/MIM_GP.cpp"

surrogate_parallel <- function(X_ls, X_type, l, fun, kernel, test_X) {
  ans_list <- progressr::with_progress({
    p <- progressr::progressor(along = X_ls)

    future.apply::future_lapply(1:length(X_ls), function(i) {
      worker_env <- environment(surrogate)
        if(!exists(".mim_gp_loaded", envir = worker_env, inherits = F)) {
          Rcpp::sourceCpp(cpp_file, env = worker_env, rebuild = F, showOutput = F, verbose = F)
          assign(".mim_gp_loaded", TRUE, envir = worker_env)
        }

      timing <- system.time({
        ans <- surrogate(X = X_ls[[i]], X_type = X_type, l = l, fun = fun, seed = i, kernel = kernel, test_X = test_X)
      }) 
      p()
      list(result = ans, elapsed_time = unname(timing["elapsed"]))
    }, future.seed = TRUE, future.packages = c("Rcpp", "nloptr", "rkriging", "MOFAT"))
  })

  mse_mat = matrix(vapply(ans_list, function(x) x$result$mse, numeric(1)), ncol = 1)
  nmse_mat = matrix(vapply(ans_list, function(x) x$result$nmse, numeric(1)), ncol = 1)
  norm_t_mat = do.call(rbind, lapply(ans_list, function(x) x$result$norm_t))
  elapsed_time = vapply(ans_list, function(x) x$elapsed_time, numeric(1))
  active_col_ls = lapply(ans_list, function(x) x$result$active_col)

  list(mse_mat = mse_mat,nmse_mat=nmse_mat,norm_t_mat = norm_t_mat, elapsed_time = elapsed_time,
  active_col_ls = active_col_ls)
}