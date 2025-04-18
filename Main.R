require(Matrix)
require(fda)
require(MASS)
require(SparseM)
require(robustbase)

rob.scale <- function(x, y, tuning.c = 1e-03, max.iter = 500, pen.reg = 2, pen.scale = 2, lambda = NULL,
                      toler = 1e-08, v = 2, type = "rob"){

## x is the univariate ordered predictor variable
## y is the response variable
## tuning.c is the tuning constant in the local quadratic approximation of the absolute loss about the origin
 # it is needed to estimate the regression function
## max.iter is the maximum number of iterations for the computing algorithms
## pen.reg is the order of the penalty for the estimation of the regression function
 # by default equal to 2 corresponding to cubic smoothing splines
## pen.scale is the order of the penalty of the estimation of the (log-)scale of the error
## lambda is the smoothing parameter for the estimation of the (log-)scale
 # NULL indicates a grid search procedure as in Kalogridis and Van Aelst (2024), otherwise specify exact value
## v is the degrees of freedom of the t-distribution and therefore the tuning parameter of the objective function
 # see Kalogridis and Van Aelst (2025) for details
  
  
  l1.smsp <- function(x, y, tuning = 1e-03, maxit = 500, m = 2, 
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
      fit.r <- l1.step(X = b.basis.e,  y = y, lambda = lambda, Pen.matrix = Pen.matrix,
                       tuning = tuning, resids.in = resids.in, maxit = maxit)
      GCV.scores <-  mean( fit.r$weights*((fit.r$resids)^2) /((1-mean(fit.r$hat.values))^2) )
      # GCV.scores <- scaleTau2(fit.r$resids/(1-fit.r$hat.values))
      return(GCV.scores)
    }
    cand.values <- exp(seq(log(1e-12), log(3), len = 30))
    GCV.values <- sapply(cand.values, FUN = GCV)
    lambda1 = cand.values[which.min(GCV.values)]
    
    lambda1 <- ifelse(is.null(lambda), lambda1, lambda)
    fit.rf <- l1.step(X = b.basis.e,  y = y, lambda = lambda1, Pen.matrix = Pen.matrix, 
                        tuning = tuning.c, resids.in = resids.in, maxit = max.iter)
    
    mu <- as.vector(b.basis.e%*%fit.rf$beta.hat)
    
    return(list(x = x, mu = mu, weights = fit.rf$weights, Pen.matrix = Pen.matrix, fitted = fit.rf$fitted, 
                lambda = lambda1, ic = fit.rf$ic, resids = fit.rf$resids, tol.f = fit.rf$check))
  }
  
  loss.f <- function(x) ((v+1)/2)*log(1+exp(2*x)/v)-x # general t
  psi.f <- function(x) ((v+1)/v)*exp(2*x)/(1+exp(2*x)/v)-1
  
  n <- length(y)
  x.c <- (x-min(x))/(max(x)-min(x))
  b.basis <- create.bspline.basis(rangeval = c(x.c[1], x.c[n]), breaks = unique(x.c), norder = 2*pen.scale)
  B <- eval.basis(b.basis, x.c)
  Pen.matrix <- bsplinepen(b.basis, Lfdobj= pen.scale)
  if(type == "rob"){
    reg.est <- l1.smsp(x,y)
    Y.log.trans <- log(abs(y-reg.est$mu))
  } else{
    # reg.est <- smooth.spline(x,y)
    GCV.l <- function(lambda){
      coef.est <- solve(t(B)%*%B+lambda*Pen.matrix, (t(B)%*%y))
      fitted <- c(B%*%coef.est)
      resids <- y-fitted
      hat.tr <-  diag(B%*%solve(t(B)%*%B+lambda*Pen.matrix, t(B)))
      GCV.crit <- mean(resids^2/(1-hat.tr)^2)
      return(GCV.crit)
    }
    cand.values <- exp(seq(log(1e-12), log(3), len = 30))
    GCV.values <- sapply(cand.values, FUN = GCV.l)
    lambda1 = cand.values[which.min(GCV.values)]
    
    coef.est <- solve(t(B)%*%B+lambda1*Pen.matrix, (t(B)%*%y))
    fitted <- c(B%*%coef.est)
    Y.log.trans <- log(abs(y-fitted))
  }
  
  pseudo.step <- function(B, y, tol = toler, maxiter = max.iter, lambda){
    
    mu.in <- 0
    ic = 0
    istop = 0
    
    while(istop == 0 & ic <= maxiter){
      ic = ic + 1
      pseudo.up <- mu.in + psi.f(y-mu.in)/2
      # pseudo.prelim.av <- apply(pseudo.up, 2, FUN = mean)
      M1 <- t(B)%*%B+ lambda*Pen.matrix
      M2 <-  t(B)%*%pseudo.up
      beta.up = SparseM::solve(M1, M2)
      mu.up <- B%*%beta.up
      check = max( abs(mu.up-mu.in) ) 
      if(check < toler){istop =1}
      mu.in <- mu.up
    }
    pseudo.f <- pseudo.up
    return(list(mu = mu.up, ic = ic, check = check, pseudo = pseudo.f))
  }
  
  if(is.null(lambda)){
    GCV <- function(lambda){
      fit.l <- pseudo.step(B=B, y = Y.log.trans, lambda = lambda)
      hat.tr <- diag( B%*%solve(t(B)%*%B+lambda*Pen.matrix, t(B)) )
      GCV.crit <- scaleTau2((fit.l$pseudo-fit.l$mu)/(1-hat.tr))
      return(GCV.crit)
    }
    cand.values <- exp(seq(log(1e-12), log(3), len = 30))
    GCV.values <- sapply(X=cand.values, FUN=GCV)
    cand.opt <- cand.values[which.min(GCV.values)]
    
    fit.f <-   pseudo.step(B, Y.log.trans, lambda = cand.opt)
  } else{
    fit.f <- pseudo.step(B, Y.log.trans, lambda = lambda)
  }
  f.int <- function(x) integrate(function(s) s^2/(x^2*v+s^2)*dnorm(s), -Inf, Inf)$value - 1/(v+1)
  cor.factor <- uniroot(f.int, interval = c(0, 1.5))$root
  mu.f <- exp(fit.f$mu)/cor.factor
  if(type == "rob"){
    reg.est = reg.est$mu
  } else{
    # reg.est = predict(reg.est, x)$y
    reg.est <- c(B%*%coef.est)
  }
  resids <- y - reg.est
  
  return(list(scale.est = mu.f, reg.est = reg.est, resids= resids, st.resids = resids/mu.f,
              lambda.scale = cand.opt,ic = fit.f$ic, check = fit.f$check))
}

