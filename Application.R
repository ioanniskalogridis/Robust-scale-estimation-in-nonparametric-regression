require(saqgetr)
require(tidyverse)
require(worldmet)
require(lubridate)
require(openair)

site.codes <- importMeta(source = "aurn", all = FALSE, year = NA, duplicate = FALSE)
kc1 <- importAURN(site = "SCN2", year = 2024) #nott OK, 
data <- select(kc1, nox, air_temp)
data <- data[complete.cases(data),]
data <- data[order(data$air_temp),]
dim(data)
y <- data$nox
x <- data$air_temp

par(mar = c(6,6.5,2,2), mgp = c(4.5, 1.5, 0))
plot(x,y, cex.axis = 2.5, cex.lab = 2.5, cex = 2.5, pch = 19, col = "gray",
     ylab = "NOx (PPB)", xlab = "Temperature (\u00B0C)") ; grid()
fit.ls1 <- rob.scale(x, y, v = 200, type = "ls")
fit.rob1 <- rob.scale(x, y, v = 1, type = "rob")

plot(x, fit.ls1$scale.est, lwd = 3, col = "black", type = "l", ylim = c(0, 50),
     lty = 2,
     cex.axis = 2.5, cex.lab = 2.5, cex = 2.5, ylab = "Estimates", xlab = "Temperature (\u00B0C)")
lines(x, fit.rob1$scale.est, lwd = 3, col = "blue"); grid()

cut.off <- 2.5
sum(abs(fit.ls1$st.resids)>cut.off)
sum(abs(fit.rob1$st.resids)>cut.off)
plot(x,y, cex.axis = 2.5, cex.lab = 2.5, cex = 2.5, pch = 19, col = "gray",
     ylab = "NOx (PPB)", xlab = "Temperature (\u00B0C)") ; grid()
points(x[abs(fit.rob1$st.resids)>cut.off], y[abs(fit.rob1$st.resids)>cut.off], 
       col = "red", cex = 2.5, pch = 19)

x.r <- x[which(fit.rob1$st.resids<=cut.off)]
y.r <- y[which(fit.rob1$st.resids<=cut.off)]
plot(x.r,y.r, ylim = c(0, 100))
fit.ls2 <- rob.scale(x.r, y.r, v = 200, type = "ls")
fit.rob2 <- rob.scale(x.r, y.r, v = 1, type = "rob")

plot(x.r, fit.ls2$scale.est, lwd = 3, col = "black", type = "l", ylim = c(0, 50),
     lty = 2,
     cex.axis = 2.5, cex.lab = 2.5, cex = 2.5, ylab = "Estimates", xlab = "Temperature (\u00B0C)")
lines(x.r, fit.rob2$scale.est, lwd = 3, col = "blue"); grid()

