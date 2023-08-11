library(CVXR)
max_smooth <- function(x) {
    return(0.5 * (sqrt(x * x + 1e-20) + x))
}

# Shannon Entropy
sentropy <- function(w) {
    if (any(w) < 0) {
        idx <- which(w < 0)
        if (all(abs(w[idx]) < 1e-04))
            w[idx] <- abs(w[idx]) else stop("\nShannon Entropy only valid for positive weights")
    }
    return(-sum(w * log(w)))
}

# Cross Entropy Approximation of Golan, Judge and Miller (1996)
# w = reference weights
# q = portfolio weights
centropy <- function(w, q) {
    if (any(q) < 0) {
        idx <- which(q < 0)
        if (all(abs(q[idx]) < 1e-04))
            q[idx] <- abs(q[idx]) else stop("\nCross Entropy weights q must be positive")
    }
    return(sum((1/q) * (w - q)^2))
}


fun_mad <- function(w, x) {
    d <- x %*% w
    d <- scale(d, scale = FALSE)
    return(mean(abs(d)))
}

fun_lpm <- function(w, x, threshold, moment) {
    if (threshold == 999) {
        x <- scale(x, scale = FALSE)
        threshold <- 0
    }
    return((mean(pmax(x %*% w - threshold, 0)^moment))^(1/moment))
}

fun_upm <- function(w, x, threshold, moment) {
    if (threshold == 999) {
        x <- scale(x, scale = FALSE)
        threshold <- 0
    }
    return((mean(pmax(threshold - x %*% w, 0)^moment))^(1/moment))
}

fun_minmax <- function(w, x) {
    return(min(x %*% w))
}

fun_var <- function(w, x) {
    d <- (x %*% w)
    d <- scale(d, scale = FALSE)
    return(mean(d^2))
}

fun_cvar <- function(w, x, VaR = NULL, alpha = 0.05) {
    if (is.null(VaR)) {
        n <- dim(x)[1]
        Rp <- (x %*% w)
        sorted <- sort(Rp)
        n.alpha <- floor(n * alpha)
        VaR <- sorted[n.alpha]
    }
    d <- VaR - (x %*% w)
    f <- -VaR + mean(max_smooth(d))/alpha
    return(f)
}

fun_cumret <- function(w, x) {
    f <- cumsum(as.numeric(x %*% w))
    return(f)
}

fun_drawdown <- function(w, x, lower = TRUE) {
    x <- fun_cumret(w, x)
    n <- length(x)
    if (lower) {
        f <- sapply(1:n, FUN = function(i) max(x[1:i]) - x[i])
    } else {
        f <- sapply(1:n, FUN = function(i) x[i] - min(x[1:i]))
    }
    return(f)
}

fun_cdar <- function(w, x, DaR = NULL, alpha = 0.05, lower = TRUE) {
    port <- x %*% w
    if (!lower) port <- -port
    dd <- fun_drawdown(w, x, lower)
    n <- length(dd)
    probability <- rep(1/n, n)
    xdd <- sort.int(dd, index.return = TRUE)
    sidx <- xdd$ix
    # key step: reverse so that we work with 0:alpha (left tail)
    xdd <- rev(xdd$x)
    cp <- cumsum(probability[sidx[1:n]])
    intalpha <- max(max(which(cp <= alpha)), 1)
    if (is.null(DaR)) DaR <- as.numeric(xdd[intalpha])
    intalpha2 <- max(max(which(cp < alpha)), 1)
    intalpha3 <- max(min(which(cp >= alpha)), 1)
    CDaRplus <- sum(probability[sidx[1:intalpha2]]/sum(probability[sidx[1:intalpha2]]) * xdd[1:intalpha2])
    lambda <- (sum(probability[1:intalpha3]) - alpha)/(alpha)
    f <- as.numeric(lambda * DaR + (1 - lambda) * CDaRplus)
    return(f)
}

riskfun <- function(weights, x, risk = c("mad", "ev", "minimax", "es", "cdar", "lpm"), alpha = 0.05,
                    moment = 1, threshold = 0, VaR = NULL, DaR = NULL, lower = TRUE) {
    xrisk <- match.arg(tolower(risk[1]), c("mad", "ev", "minimax", "es", "cdar", "lpm"))
    ans <- switch(tolower(xrisk),
                  mad = fun_mad(weights, x),
                  ev = fun_var(weights, x),
                  minimax = fun_minmax(weights, x),
                  es = fun_cvar(weights, x, VaR = VaR, alpha = alpha),
                  lpm = fun_lpm(weights, x, moment = moment, threshold = threshold),
                  cdar = fun_cdar(weights, x, DaR = DaR, alpha = alpha, lower = lower))
    return(abs(ans))
}
