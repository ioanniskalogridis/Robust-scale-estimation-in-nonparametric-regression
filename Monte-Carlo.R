require(fda)
require(sn)
require(EnvStats)

nrep <- 1000
n <- 120
m.scale.ls <- matrix(NA, nrow = n, ncol = nrep)
m.scale.rob <- matrix(NA, nrow = n, ncol = nrep)

mse.scale.ls <- rep(NA, nrep)
mse.scale.rob <- rep(NA, nrep)

for(k in 1:nrep){
  print(k)
  y <- rep(0,n)
  x <- 1:n/n
  for(i in 1:n){
    y[i] <- cos(2*pi*x[i])+ 0.5*exp(x[i])*rnorm(1)
  }
  # r.samp <- sample(1:n, floor(0.05*n)) # 0.05 and 0.1 contamination levels
  # y[r.samp] <- 5*y[r.samp]
  fit.rob <- rob.scale(x, y, v = 2)
  fit.ls <- rob.scale(x, y, v = 50)
  
  mse.scale.ls[k] <- mean((fit.ls$mu-0.5*exp(x))^2)
  mse.scale.rob[k] <- mean((fit.rob$mu-0.5*exp(x))^2)
  
  m.scale.ls[, k] <- fit.ls$mu
  m.scale.rob[, k] <- fit.rob$mu
}
mean(mse.scale.ls, na.rm = TRUE);mean(mse.scale.rob, na.rm = TRUE)

matplot(x,m.scale.ls, type = "l", lwd = 3, lty = 1, cex.axis = 2.5, cex.lab = 2.5, cex = 3, ylab = "y",
        col = "gray") ; grid()
curve(0.5*exp(x), add= TRUE, lwd = 3, col = "black")
matplot(x,m.scale.rob, type = "l", lwd = 3, lty = 1, cex.axis = 2.5, cex.lab = 2.5, cex = 3, ylab ="y",
        col = "gray", ylim = c(0.3, 2.5)); grid()
curve(0.5*exp(x), add= TRUE, lwd = 3)

