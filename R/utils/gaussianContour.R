gaussianContour <- function(mu, Sigma, val) {
    # assumes 2D
    stopifnot(all(dim(Sigma)==c(2,2)))
    
    npoints <- 100
#     A <- log(det(Sigma*2*pi))
#     tmp <- solve(Sigma, mu)
#     const <- sum(mu*tmp) + A - log(val)
    const <- val
    
    circ <- cbind(sin(seq(0,2*pi,len=npoints)),cos(seq(0,2*pi,len=npoints)))*const

    R <- chol(solve(Sigma))
    
    out <- t(solve(R, t(circ))) +  matrix(mu, npoints, 2, byrow = TRUE)

    return(out)
}