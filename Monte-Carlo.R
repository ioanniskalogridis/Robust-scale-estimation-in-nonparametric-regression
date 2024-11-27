require(fda)
require(sn)
require(EnvStats)

nrep <- 500
n <- 100
m.scale.ls <- matrix(NA, nrow = n, ncol = nrep)
m.scale.rob1 = m.scale.rob2 = m.scale.rob5 = m.scale.rob10 <- matrix(NA, nrow = n, ncol = nrep)

mse.scale.ls <- rep(NA, nrep)
mse.scale.rob1 = mse.scale.rob2 = mse.scale.rob5 = mse.scale.rob10 <- rep(NA, nrep)

f1 <- function(x) 0.5*exp(x)

for(k in 1:nrep){
  print(k)
  y <- rep(0,n)
  x <- 1:n/n
  y <- cos(2*pi*x) + f1(x)*rnorm(n)
  # y <- cos(4*pi*x) + f1(x)*rnorm(n)
  # par(mar = c(6,6.5,2,2), mgp = c(4.5, 1.5, 0))
  # plot(x, y, cex = 2.5, col = "gray", pch = 19, ylab = "Y", xlab = "x", cex.lab = 2.5,
       # cex.axis = 2.5) ; grid()
  # curve(cos(4*pi*x), add = TRUE, lwd = 4, col = "black")
  # r.samp <- sample(1:n, floor(0.05*n)) # 0.05, 0.1 and 0.2 contamination levels
  # y[r.samp] <- 5*y[r.samp]
  fit.rob1 <- rob.scale(x, y, v = 1)
  fit.rob2 <- rob.scale(x, y, v = 2)
  fit.rob5 <- rob.scale(x, y, v = 5)
  fit.rob10 <- rob.scale(x, y, v = 10)
  fit.ls <- rob.scale(x, y, v = 50)
  
  mse.scale.ls[k] <- mean((fit.ls$mu-f1(x))^2)
  mse.scale.rob1[k] <- mean((fit.rob1$mu-f1(x))^2)
  mse.scale.rob2[k] <- mean((fit.rob2$mu-f1(x))^2)
  mse.scale.rob5[k] <- mean((fit.rob5$mu-f1(x))^2)
  mse.scale.rob10[k] <- mean((fit.rob10$mu-f1(x))^2)
  
  m.scale.ls[, k] <- fit.ls$mu
  m.scale.rob1[, k] <- fit.rob1$mu
  m.scale.rob2[, k] <- fit.rob2$mu
  m.scale.rob5[, k] <- fit.rob5$mu
  m.scale.rob10[, k] <- fit.rob10$mu
}

mean(mse.scale.rob1, na.rm = TRUE)*100;sd(mse.scale.rob1, na.rm = TRUE)/sqrt(nrep)*100
mean(mse.scale.rob2, na.rm = TRUE)*100; sd(mse.scale.rob2, na.rm = TRUE)/sqrt(nrep)*100
mean(mse.scale.rob5, na.rm = TRUE)*100; sd(mse.scale.rob5, na.rm = TRUE)/sqrt(nrep)*100
mean(mse.scale.rob10, na.rm = TRUE)*100; sd(mse.scale.rob10, na.rm = TRUE)/sqrt(nrep)*100
mean(mse.scale.ls, na.rm = TRUE)*100;  sd(mse.scale.ls, na.rm = TRUE)/sqrt(nrep)*100

matplot(x, m.scale.ls, type = "l", lwd = 3, lty = 1, cex.axis = 2.5, cex.lab = 2.5, cex = 3, ylab = "y",
        col = "gray", ylim = c(0.3, 2.5)) ; grid()
curve(f1(x), add= TRUE, lwd = 3, col = "black")
matplot(x, m.scale.rob2, type = "l", lwd = 3, lty = 1, cex.axis = 2.5, cex.lab = 2.5, cex = 3, ylab ="y",
        col = "gray", ylim = c(0.3, 2.5)); grid()
curve(f1(x), add= TRUE, lwd = 3, col = "black")

