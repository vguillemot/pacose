ridge.cv <- function (X, y, lambda = NULL, scale = TRUE, k = 10, plot.it = FALSE, nlambda = 100) 
{
    if (is.null(lambda) == TRUE) {
        ss <- seq(-10, -1, length = nlambda)
        ss <- 10^ss
        n <- nrow(X)
        nn <- n - floor(n/k)
        lambda <- ss * nn * ncol(X)
    }
    cv <- rep(0, length(lambda))
    n <- nrow(X)
    all.folds <- split(sample(seq_len(n)), rep(seq_len(k), length = n))
    for (i in seq_len(k)) {
        omit <- all.folds[[i]]
        Xtrain = X[-omit, , drop = FALSE]
        ytrain = y[-omit]
        Xtest = X[omit, , drop = FALSE]
        ytest = y[omit]
        ll <- lm.ridge(ytrain ~ Xtrain, scale = scale, lambda = lambda)
        coef.ll <- coef(ll)
        pred <- t(matrix(coef.ll[, 1], nrow = length(lambda), 
            ncol = length(ytest))) + Xtest %*% t(coef.ll[, -1])
        cv <- cv + colSums((pred - ytest)^2)
    }
    cv <- cv/n
    lambda.opt <- lambda[which.min(cv)]
    if (plot.it == TRUE) {
        plot(lambda, cv, type = "l")
    }
    rr <- lm.ridge(y ~ X, scale = scale, lambda = lambda.opt)
    coefficients <- coef(rr)
    intercept <- coefficients[1]
    coefficients <- coefficients[-1]
    return(list(intercept = intercept, coefficients = coefficients, 
        lambda.opt = lambda.opt))
}

