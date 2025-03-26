ridge.net <-
function (X, lambda = NULL, plot.it = FALSE, scale = TRUE, k = 10,
    verbose = FALSE, cv.method="HKB")
{
    if (is.null(lambda) == TRUE) {
        ss <- seq(-10, -1, length = 1000)
        ss <- 10^ss
        n <- nrow(X)
        nn <- n - floor(n/k)
        lambda <- ss * nn * ncol(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    X <- scale(X, scale = scale)
    B <- matrix(0, nrow = p, ncol = p)
    lambda.opt <- rep(0, p)
    if (verbose) cat(paste("Performing local ridge regressions\n"))
    if (verbose) cat(paste("Vertex no "))
    for (i in 1:p) {
        if ((i/10) == floor(i/10)) {
            if (verbose) cat(paste(i, "..."))
        }
        noti <- (1:p)[-i]
        yi <- X[, i]
        Xi <- X[, noti]
        if (cv.method == "CV" ) lambda.opt.i <- ridge.cv(Xi, yi, lambda = lambda, scale = scale, plot.it = plot.it, k = k)$lambda.opt
        if (cv.method == "HKB") lambda.opt.i <- lm.ridge(yi~Xi, lambda = 0, scale = scale)$kHKB
        rr <- lm.ridge(yi ~ Xi, scale = scale, lambda = lambda.opt.i)
        B[i, -i] <- coef(rr)[-1]
        lambda.opt[i] <- lambda.opt.i
    }
    if (verbose) cat("\n")
    pcor <- Beta2parcor(B, verbose = verbose)
    return(list(pcor = pcor,lambda.opt=lambda.opt))
}

