pos <- c(4,2)
transition <- matrix(c(1,1,0,1),2,2)
Q <- diag(2)*0.3

emission <- matrix(c(1,0,0,1),2,2)
R <- matrix(c(1,0.5,0.5,1),2,2)*0.03

nobs    <- 20
y       <- matrix(0, nobs, 2)
x       <- matrix(0, nobs, 2)

for(ii in 1:nobs) {
    pos <- transition %*% pos + t(mvtnorm::rmvnorm(1,sigma=Q))
    x[ii,] <- pos
    y[ii,] <- emission %*% pos + t(mvtnorm::rmvnorm(1,sigma=R))
}

tmp <- kalmanFilter(y, matrix(as.vector(emission),nobs,4,byrow=TRUE), transition, R, Q)

# -------do sinusoidal system----------

transition <- matrix(c(1,0,0,1),2,2)
Q <- diag(2)*0.00003

emission <- matrix(c(1,0,0,1),2,2)
R <- 0.2

nobs    <- 200
x       <- sin((1:nobs)*2*pi/80)
y       <- x + rnorm(200,sd =R^2)
H       <- cbind(1, 1:nobs)

tmp <- kalmanFilter(matrix(y, ncol=1), H, transition, matrix(R,1,1), Q, bPlot=FALSE)
theta <- tmp$x
yhat <- rowSums(H * theta)
opar <- par(mfrow = c(2,1))
plot(y)
lines(x, col = "red")
lines(yhat, col = "steelblue", lty=3)
plot(1:nobs, theta[,2])
