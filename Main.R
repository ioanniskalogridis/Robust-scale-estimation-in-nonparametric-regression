require(Matrix)
require(fda)
require(MASS)
require(SparseM)

rob.scale <- function(x, y, tuning.c = 1e-04, max.iter = 500, pen.reg = 2, pen.scale = 2, lambda = NULL,
                      toler = 1e-08){
  
  l1.smsp <- function(x, y, tuning = 1e-04, maxit = 500, m = 2, 
                        lambda = NULL){
    
    rho.abs <-  function(x, tuning = 1e-04) ifelse(abs(x) <= tuning, 0.5*x^2/tuning, abs(x))
    psi.abs <- function(x, tuning = 1e-04) ifelse(abs(x) <= tuning, x/tuning, sign(x))
    weight.abs <- function(x, tuning = 1e-04) ifelse(abs(x) <= tuning, 1/tuning, sign(x)/x)
    
    l1.step <- function(X, y, tuning = tuning.c , tol = toler, lambda, Pen.matrix, resids.in, maxit = max.iter){
      
      ic = 0
      istop = 0
      
      while(istop == 0 & ic <= maxit){
        ic = ic + 1
        weights.prelim <- weight.abs(resids.in, tuning = tuning)
        X.s <- as.matrix(scale(t(X), center = FALSE, scale = 1/weights.prelim))
        M1 <-  X.s%*%X + 2*lambda*Pen.matrix
        M2 <- X.s%*%as.matrix(y)
        v1 = SparseM::solve( M1, M2)
        resids1 <- as.vector(y-X%*%v1)
        check = max( abs(resids1-resids.in)) 
        if(check < tol){istop =1}
        resids.in <- resids1
      }
      weights1 = as.vector(weight.abs(resids1, tuning = tuning) )
      fitted.values = as.vector(X%*%v1)
      resids = y - fitted.values
      hat.values = diag(X%*%SparseM::solve(M1, X.s))
      # hat.values = sum( t(SparseM::solve(M1, X.s)) * X )/length(y)
      return(list(resids = resids, beta.hat = v1, hat.values = hat.values, 
                  weights = weights1, fitted = fitted.values, ic = ic, check = check ) )
    }
    n <- length(y)
    
    x.t <- (x-min(x))/(max(x)-min(x))
    
    b.basis <- create.bspline.basis(rangeval = c(0, 1), breaks = unique(x.t), norder = 2*pen.reg)
    b.basis.e <- eval.basis(b.basis, x.t)
    
    if(pen.reg==1){
      P <- diff(diag(n), differences = 1)
      Pen.matrix <- t(P)%*%diag(1/(diff(x, differences = 1)))%*%P
    } else {
      Pen.matrix <- bsplinepen(b.basis, Lfdobj = pen.reg) }
    
    fit.in <- smooth.spline(x, y)
    resids.in <- y - predict(fit.in, x.t)$y
    GCV <- function(lambda){
      fit.r <- quan.step(X = b.basis.e,  y = y, lambda = lambda, Pen.matrix = Pen.matrix,
                         alpha = alpha, tuning = tuning, resids.in = resids.in, maxit = maxit)
      GCV.scores <-  mean( fit.r$weights*((fit.r$resids)^2) /((1-fit.r$hat.values)^2) )
      return(GCV.scores)
    }
    cand.values <- exp(seq(log(1e-12), log(3), len = 50))
    GCV.values <- sapply(cand.values, FUN = GCV)
    lambda1 = cand.values[which.min(GCV.values)]
    
    lambda1 <- ifelse(is.null(lambda), lambda1, lambda)
    fit.rf <- quan.step(X = b.basis.e,  y = y, lambda = lambda1, Pen.matrix = Pen.matrix, 
                        tuning = tuning.c, resids.in = resids.in, maxit = max.iter)
    
    mu <- as.vector(b.basis.e%*%fit.rf$beta.hat)
    
    return(list(x = x, mu = mu, weights = fit.rf$weights, Pen.matrix = Pen.matrix, fitted = fit.rf$fitted, 
                lambda = lambda1, ic = fit.rf$ic, resids = fit.rf$resids, tol.f = fit.rf$check))
  }
  
  
  
}


quan.smsp <- function(x, y, tuning = 1e-04, maxit = 500, m = 2, 
                      lambda = NULL){
  
  rho.abs <-  function(x, tuning = 1e-04) ifelse(abs(x) <= tuning, 0.5*x^2/tuning, abs(x))
  psi.abs <- function(x, tuning = 1e-04) ifelse(abs(x) <= tuning, x/tuning, sign(x))
  weight.abs <- function(x, tuning = 1e-04) ifelse(abs(x) <= tuning, 1/tuning, sign(x)/x)
  
  quan.step <- function(X, y, tuning, alpha,  tol = 1e-08, lambda, Pen.matrix, resids.in, maxit = 50){
    
    ic = 0
    istop = 0
    
    while(istop == 0 & ic <= maxit){
      ic = ic + 1
      weights.prelim <- weights.quan(resids.in, alpha = alpha, tuning = tuning)
      X.s <- as.matrix(scale(t(X), center = FALSE, scale = 1/weights.prelim))
      M1 <-  X.s%*%X + 2*lambda*Pen.matrix
      M2 <- X.s%*%as.matrix(y)
      v1 = SparseM::solve( M1, M2)
      # v1 = ginv(M1)%*%M2
      resids1 <- as.vector(y-X%*%v1)
      check = max( abs(resids1-resids.in)) 
      if(check < tol){istop =1}
      resids.in <- resids1
    }
    weights1 = as.vector(weights.quan(resids1, alpha = alpha, tuning = tuning) )
    fitted.values = as.vector(X%*%v1)
    resids = y - fitted.values
    hat.values = diag(X%*%SparseM::solve(M1, X.s))
    # hat.values = sum( t(SparseM::solve(M1, X.s)) * X )/length(y)
    
    return(list(resids = resids, beta.hat = v1, hat.values = hat.values, 
                weights = weights1, fitted = fitted.values, ic = ic, check = check ) )
  }
  
  n <- length(y)
  
  x.t <- (x-min(x))/(max(x)-min(x))
  
  b.basis <- create.bspline.basis(rangeval = c(0, 1), breaks = unique(x.t), norder = 2*m)
  b.basis.e <- eval.basis(b.basis, x.t)
  
  if(m==1){
    P <- diff(diag(n), differences = 1)
    Pen.matrix <- t(P)%*%diag(1/(diff(x, differences = 1)))%*%P
  } else {
    Pen.matrix <- bsplinepen(b.basis, Lfdobj = m) }
  
  fit.in <- smooth.spline(x, y)
  resids.in <- y - predict(fit.in, x.t)$y
  GCV <- function(lambda){
    fit.r <- quan.step(X = b.basis.e,  y = y, lambda = lambda, Pen.matrix = Pen.matrix,
                       alpha = alpha, tuning = tuning, resids.in = resids.in, maxit = maxit)
    GCV.scores <-  mean( fit.r$weights*((fit.r$resids)^2) /((1-fit.r$hat.values)^2) )
    return(GCV.scores)
  }
  cand.values <- exp(seq(log(1e-12), log(3), len = 50))
  GCV.values <- sapply(cand.values, FUN = GCV)
  lambda1 = cand.values[which.min(GCV.values)]
  
  lambda1 <- ifelse(is.null(lambda), lambda1, lambda)
  fit.rf <- quan.step(X = b.basis.e,  y = y, lambda = lambda1, Pen.matrix = Pen.matrix,
                      alpha = alpha, tuning = tuning, resids.in = resids.in, maxit = maxit)
  
  mu <- as.vector(b.basis.e%*%fit.rf$beta.hat)
  
  return(list(x = x, mu = mu, weights = fit.rf$weights, Pen.matrix = Pen.matrix, fitted = fit.rf$fitted, 
              lambda = lambda1, ic = fit.rf$ic, resids = fit.rf$resids, tol.f = fit.rf$check))
}