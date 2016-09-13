kalmanFilter <- function(y, H, A, R, Q, prior=NULL, bPlot=TRUE, plotlim=NULL) {
    
    # y   - (n x d) matrix of input observations
    # H   - (n x (k * d)) matrix of likelihood multipliers. Unravelled col-wise
    # A   - (k x k) matrix of mean state transition.
    # R   - (d x d) matrix for likelihood covariance
    # Q   - (k x k) matrix for state transition covariance
    # prior - list of (mean, var) for prior of x0
    
    # no care has been taken to use a numerically stable version of the algorithm
    
    assert(all(sapply(list(y, H, A, R, Q), is.matrix)), "inputs y, H, A, R, Q must all be matrices")
    
    state   <- list()
    obs     <- list()
    n       <- dim(y)[1]
    d       <- dim(y)[2]
    k       <- dim(A)[1]
    X       <- matrix(0, n, k)
    allm    <- matrix(0, n, k)
    allP    <- matrix(0, n, k*k)
    
    assert(all(dim(H)==c(n, d*k)), "H must be of size (n, (k*d)) = (%d, %d)", n, d*k)
    assert(all(dim(A)==c(k,k)), "A must be of size (k, k) = (%d, %d)", k, k)
    assert(all(dim(R)==c(d,d)), "R must be of size (d, d) = (%d, %d)", d, d)
    assert(all(dim(Q)==c(k,k)), "Q must be of size (k, k) = (%d, %d)", k, k)
    
    if(is.null(prior)) {
        prior <- list(mean=rep(0, k), var=diag(k)*1e3)
    }
    assert(all(names(prior) %in% c("mean", "var")), "prior must have mean and var fields")

    
    m <- list()
    P <- list()
    loop <- list()
    m$minus <- A %*% prior$mean
    P$minus <- A %*% prior$var %*% t(A) + Q
    
    if(bPlot && is.null(plotlim)) {
        plotlim <- list()
        plotlim$rngx <- diff(range(y[,1]))
        plotlim$rngy <- diff(range(y[,2]))
        plotlim$x <- c(min(y[,1])-0.1*plotlim$rngx, max(y[,1]) + 0.1*plotlim$rngx)
        plotlim$y <- c(min(y[,2])-0.1*plotlim$rngy, max(y[,2]) + 0.1*plotlim$rngy)
    }
    
    for(ii in 1:n) {
        loop$H <- matrix(H[ii,], d, k)
        loop$v <- y[ii,] - loop$H %*% m$minus
        loop$S <- loop$H %*% P$minus %*% t(loop$H) + R
        loop$K <- P$minus %*% t(loop$H) %*% solve(loop$S)
        m$curr <- m$minus + loop$K %*% loop$v
        P$curr <- P$minus - loop$K %*% loop$S %*% t(loop$K)
        
        if(bPlot && ii > 1) {
            plot(y[1:(ii-1),1],y[1:(ii-1),2], xlim=plotlim$x, ylim=plotlim$y)
            lines(x[1:(ii-1),1],x[1:(ii-1),2])
            lines(gaussianContour(m$minus, P$minus, 2), lty=2, col="red")
            Sys.sleep(time = 0.3)
            points(y[ii,1], y[ii,2], pch=3, cex=3, col="steelblue")
            lines(gaussianContour(m$curr, P$curr, 2), lty=3)
            cat(sprintf("Press [enter] to continue..."))
            readline()
        }
        
        X[ii,]      <- m$curr
        allP[ii,]   <- as.vector(P$curr)
        
        m$minus <- A %*% m$curr
        P$minus <- A %*% P$curr %*% t(A) + Q
    }
    return(list(x=X, P = allP))
}
            
            
