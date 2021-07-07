
# Comment these out, if Rcpp/C++ compiler are not present.
library(Rcpp)
sourceCpp(file = "uniprod.cpp")

# setup rng --------------------------------
# Standard R random setup. rgamma uses this always.
# RNG <- runif
# SETSEED <- set.seed

# This needs C++ functions/compiler.
# Faster random number generator.
RNG <- prodrunif
SETSEED <- setseed
# --------------------------------------------

# prod.tozero()
prod.tozero <- function() {
    prod <- 1
    n <- 0
    while (T) {
        n <- n + 1
        prod <- prod * runif(1)
        if (prod == 0) break
    }
    n
}
# Seeding for C++ RNG
prod.set.seed <- function(seed) {
    invisible(setseed(seed))
}

# 128-bit arithmetic
# prod.mean128(samplesize=1e5, N = 200, seed = 0)
prod.mean128 <- function(samplesize = 1e6, N = 200, seed = 0) {
    invisible(prodmean128(samplesize, N, seed))
}

# 80-bit extended arithmetic
# prod.mean80(samplesize=1e6, N = 200, gamma = F, seed = 0)
prod.mean80 <- function(samplesize = 1e6, N = 200, gamma = F, seed = 0) {
    if (gamma && seed > 0) set.seed(seed) # for rgamma in C++ funcs.
    invisible(prodmean80(samplesize, N, gamma, seed))
}

# prod.mean64(samplesize=1e6, N = 200, gamma = F, seed = 0)
prod.mean64 <- function(samplesize = 1e6, N = 200, gamma = F, seed = 0) {
    if (gamma && seed > 0) set.seed(seed)
    invisible(prodmean64(samplesize, N, gamma, seed))
}

# prod.mean generates samplesize U(0, 1) products of length N and
# calculates a summary.
# prod.mean(samplesize=1e6, N=200, gamma=F, seed=1)
prod.mean <- function(samplesize = 1e6, N = 200, gamma = F, seed = 0) {
    meanprod <- 0
    if (gamma && seed > 0) set.seed(seed) # rgamma
    if (!gamma && seed > 0) SETSEED(seed)

    products <- numeric(samplesize)
    for (i in 1:samplesize) {
        if (gamma) {
            p <- exp(-rgamma(1, N))
        } else {
            p <- prod(RNG(N))
        }
        # if (N == 1) p <- p * p * p * p * p
        products[i] <- p
    }
    meanprod <- mean(products)

    w <- writeLines
    s <- sprintf
    w(s("sample size         %1.0e", samplesize))
    w(s("N                   %d", N))

    w(s("mean sample         %1.1e  %1.15e", meanprod, meanprod))
    w(s("mean distribution   %1.1e  2^-N", 2^-N))
    w(s("e^-N                %1.1e", exp(-N)))
    w(s("median sample       %1.3e", median(products)))
    w(s("median distribution %1.3e", exp(-qgamma(0.5, N))))
    w(s(
        "SE mean sample      %1.1e  by sample mean & var)",
        sd(products) / sqrt(samplesize)
    ))
    w(s(
        "SE mean sample      %1.1e  by distribution mean & var",
        sqrt((1 / 12 + 0.25)^N - 0.25^N) / sqrt(samplesize)
    ))
    w(s("log mean            %1.1f", log(meanprod)))
    w(s("mean log(prod)      %1.1f", mean(log(products))))
    w(s("median log(prod)    %1.1f", median(log(products))))
    w(s("mean geom sample    %1.5f", mean(products^(1 / N))))
    w(s(
        "(1 + 1/N)^-N        %1.5f  mean geom distribution",
        (1 + 1 / N)^-N
    ))
    w(s("1/e                 %1.5f", 1 / exp(1)))
}

# https://en.wikipedia.org/wiki/Poisson_distribution#Confidence_interval
prod.clim <- function(samplesize = 1e6, N = 200, gamma = F, seed = 0) {
    within <- 0
    lowlim <- exp(-qgamma(1 - 0.05 / 2, N))
    uplim <- exp(-qgamma(0.05 / 2, N))
    center <- exp(-qgamma(0.5, N))

    # me <- N #- 1 / 3 + (8 / 405) / N #+ 184 / (25515 * N * N)
    # lowlim <- exp(-me - 2 * sqrt(N))
    # uplim <- exp(-me + 2 * sqrt(N))
    # center <- exp(-me)

    if (gamma && seed > 0) set.seed(seed) # rgamma
    if (!gamma && seed > 0) SETSEED(seed)

    for (i in 1:samplesize) {
        if (gamma) {
            p <- exp(-rgamma(1, N))
        } else {
            p <- prod(RNG(N))
        }
        if (p >= lowlim && p <= uplim) within <- within + 1
    }
    me <- N - 1 / 3 + (8 / 405) / N #+ 184 / (25515 * N * N)

    w <- writeLines
    s <- sprintf
    w(s("sample size      %1.0e", samplesize))
    w(s("N                %d", N))
    w(s("e^-N             %1.1e", exp(-N)))
    w(s("CL low           %1.1e", lowlim))
    w(s("median           %1.4e", center))
    w(s("median est.      %1.4e", exp(-me)))
    w(s("CL high          %1.1e", uplim))
    w(s("2^-N             %1.1e", 2^-N))
    w(s("within CI        %1.4f", within / samplesize))
}

# prod.long(N) multiplies N U(0, 1) random variables and
# compares the log(result) to log(e^-N) = N.
#
# prod.long(N=1e7)
prod.long <- function(N = 1e7) {
    prod <- 1
    k <- 0
    exp2 <- 0
    k <- 0
    while (k <= N - 100) {
        k <- k + 100
        prod <- prod * prod(RNG(100))
        if (prod < 2^-800) {
            prod <- prod * 2^800
            exp2 <- exp2 + 800
        }
    }
    lprod <- log(prod) - exp2 * log(2)
    diff <- (k + lprod) / k
    w <- writeLines
    s <- sprintf
    w(s("N = %1.0f", k))
    w(s("log(prod) = -log(%1.2e x 2^-%d) = %1.1f", prod, exp2, -lprod))
    w(s("Relative difference to N = -log(e^-N) = %1.5f", -diff))
}

prod.mulToLim <- function(lim, nmax = 1000, log = F) {
    rn <- RNG(nmax)
    if (log) lim <- exp(lim)
    prod <- 1
    for (n in 1:nmax) {
        prod <- prod * rn[n]
        if (prod < lim) {
            return(n)
        }
    }
    n
}

prod.sample <- function(lim, samplesize = 10000, nmax = 1000, log = F) {
    cnt <- numeric(nmax)
    for (i in 1:samplesize) {
        k <- prod.mulToLim(lim, nmax, log)
        cnt[k] <- cnt[k] + 1
    }
    cnt
}

prod.meanLen <- function(cnt) {
    n <- length(cnt)
    sum <- sum(cnt)
    mean <- 0
    for (i in 1:n) {
        mean <- mean + i * cnt[i] / sum
    }
    mean
}

# en.wikipedia.org/wiki/Poisson_distribution#Random_drawing_from_the_Poisson_distribution
# "A simple algorithm to generate random Poisson-distributed numbers has been given by Knuth"
rpoissonKnuth <- function(lambda) {
    L <- exp(-lambda)
    k <- 0
    p <- 1
    while (TRUE) {
        k <- k + 1
        p <- p * runif(1)
        if (p < L) {
            return(k - 1)
        }
    }
}

# prod.plotPDF runs samplesize products to lim and plot the sample
#  distrubition with Poisson(-log(lim)) distribution.
#
# prod.plotPDF(lim=-745, samplesize = 1e+5, log =T)
prod.plotPDF <- function(lim, samplesize = 1e+4, log = T) {
    if (!log && lim <= 2.5e-308) {
        stop("prod.plotPDF: limit must be > 2.5e-308")
    }
    if (log && lim < -745) {
        stop("prod.plotPDF: log(limit) must be >= -745")
    }
    lambda <- -lim
    if (!log) lambda <- -log(lim)

    sample <- prod.sample(lim, samplesize, nmax = 1000, log) / samplesize
    xmax <- lambda + 3 * sqrt(lambda)
    xmin <- lambda - 3 * sqrt(lambda)
    if (xmin < 0) xmin <- 0
    pmax <- max(sample)
    plot.new()
    plot(sample,
        xlim = c(xmin, xmax), ylim = c(0, pmax),
        cex.main = 1, cex = 0.75, font.main = 1,
        xlab = paste(
            "Number of U(0,1) variables needed to reach limit = e^-",
            format(lambda, digits = 4)
        ),
        ylab = "Point  mass  probability",
        main = c(
            "Poisson/Gamma and  sample PDF distributions",
            paste("Mean length of sample product =", format(prod.meanLen(sample), digits = 4))
        )
    )
    abline(v = lambda + 1, lty = 2)
    par(new = TRUE)

    f <- function(x) dgamma(lambda, x)
    # f <- function(x) dpois(floor(x - 1), lambda)

    plot(f,
        lwd = 1.75, xlim = c(xmin, xmax),
        yaxt = "n", xaxt = "n", xlab = "", ylab = "",
        ylim = c(0, pmax)
    )
}

# prod.chisqTest tests sample and Poisson distribution as
# in plots by prod.plotPDF.
#
# prod.chisqTest(-100, samplesize = 1e+5, log = T)
prod.chisqTest <- function(lim, samplesize = 1e+5, log = F) {
    if (!log && lim <= 2.5e-308) {
        stop("prod.chisqTest: lim must be > 2.5e-308")
    }
    if (log && lim < -745) {
        stop("prod.plotPDF: log(limit) must be >= -745")
    }
    nmax <- 1000
    lambda <- -lim
    if (!log) lambda <- -log(lim)

    cnt <- prod.sample(lim, samplesize, nmax, log)
    samplemean <- prod.meanLen(cnt)
    prob <- numeric(nmax)
    for (i in 1:nmax) prob[i] <- dgamma(lambda, shape = i)
    # for (i in 1:nmax prob[i] <- dpois(i - 1, lambda)

    for (start in 1:nmax) if (prob[start] * samplesize > 25) break
    for (end in nmax:1) if (prob[end] * samplesize > 25) break

    prob[start] <- sum(prob[1:start])
    cnt[start] <- sum(cnt[1:start])
    prob[end] <- sum(prob[nmax:end])
    cnt[end] <- sum(cnt[end:nmax])
    prob <- prob[start:end]
    cnt <- cnt[start:end]
    pv <- chisq.test(x = cnt, y = NULL, p = prob)$p.value
    w <- writeLines
    s <- sprintf
    w(s("Sample mean %1.2f", samplemean))
    w(s("lambda      %1.2f", lambda))
    w(s("p-value     %1.1e", pv))
    # return(pv)
}