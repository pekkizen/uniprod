
# setup rng
library(dqrng)
RNG <- dqrunif
dqRNGkind("Xoroshiro128+")
dqset.seed(1)
# RNG <- runif
# set.seed(1)

# prod.mean samples prodcnt U(0, 1) products of length N and
# calculates a summary.
#
# prod.mean(prodcnt = 1e+6, N = 100)
prod.mean <- function(prodcnt = 1e+6, N = 100) {
    meanprod <- 0
    meanlog <- 0
    meangeom <- 0
    for (i in 1:prodcnt) {
        prod <- prod(RNG(N))
        log <- log(prod)
        # log <- -rgamma(1, N)
        # prod <- exp(log)
        meanprod <- meanprod + prod
        meanlog <- meanlog + log
        meangeom <- meangeom + prod^(1 / N)
    }
    mean <- meanprod / prodcnt
    meanlog <- meanlog / prodcnt
    meangeom <- meangeom / prodcnt

    w <- writeLines
    s <- sprintf
    w(s("rounds            %1.0e", prodcnt))
    w(s("mean(prod)        %1.2e", mean))
    w(s("log(mean(prod))   %1.1f", log(mean)))
    w(s("mean(log(prod))   %1.1f", meanlog))
    w(s("meanGeom(prod)    %1.5f", meangeom))
    w(s("(1 + 1/N)^-N      %1.5f", (1 + 1 / N)^-N))
    w(s("1/e               %1.5f", 1 / exp(1)))
    w(s("2^-N              %1.0e", 2^-N))
    w(s("2^-N/mean         %1.0e", 2^-N / mean))
    w(s("e^-N/mean         %1.0e", exp(-N) / mean))
}

# prod.lonprod(N) multiplies N U(0, 1) random variables and
# compares the log(result) to log(e^-N) = N.
#
# prod.lonprod(100000)
prod.lonprod <- function(N = 1000000) {
    prod <- 1
    k <- 0
    pow2 <- 0
    k <- 0
    while (k <= N - 100) {
        k <- k + 100
        prod <- prod * prod(RNG(100))
        if (prod < 0x1p-800) {
            prod <- prod * 0x1p+600
            pow2 <- pow2 + 600
        }
    }
    lprod <- log(prod) - pow2 * log(2)
    diff <- (k + lprod) / k
    w <- writeLines
    s <- sprintf
    w(s("N = %1.0f", k))
    w(s("log(prod) = -log(%1.2e x 2^-%d) = %1.0f", prod, pow2, -lprod))
    w(s("Relative difference to N = -log(e^-N) = %1.5f", diff))
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

prod.sample <- function(lim, rounds = 10000, nmax = 1000, log = F) {
    cnt <- numeric(nmax)
    for (i in 1:rounds) {
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

# prod.plotPDF runs prodcnt product to lim and plot the sample
#  distrubition with Poisson(-log(limit)) distribution.
#
# prod.plotPDF(lim=-745, prodcnt = 1e+5, log =T)
prod.plotPDF <- function(lim, prodcnt = 1e+4, log = T) {
    if (!log && lim <= 2.5e-308) {
        stop("prod.plotPDF: limit must be > 2.5e-308")
    }
    if (log && lim < -745) {
        stop("prod.plotPDF: log(limit) must be >= -745")
    }
    lambda <- -lim
    if (!log) lambda <- -log(lim)

    sample <- prod.sample(lim, prodcnt, nmax = 1000, log) / prodcnt
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
# prod.chisqTest(-100, rounds = 1e+5, log = T)
prod.chisqTest <- function(lim, prodcnt = 1e+5, log = F) {
    if (!log && lim <= 2.5e-308) {
        stop("prod.chisqTest: lim must be > 2.5e-308")
    }
    nmax <- 1000
    lambda <- -lim
    if (!log) lambda <- -log(lim)

    cnt <- prod.sample(lim, prodcnt, nmax, log)
    samplemean <- prod.meanLen(cnt)
    prob <- numeric(nmax)
    for (i in 1:nmax) prob[i] <- dgamma(lambda, shape = i)
    # for (i in 1:nmax prob[i] <- dpois(i - 1, lambda)

    for (start in 1:nmax) if (prob[start] * prodcnt > 25) break
    for (end in nmax:1) if (prob[end] * prodcnt > 25) break

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
