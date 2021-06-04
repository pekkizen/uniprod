
# setup rng
library(dqrng)
RNG <- dqrunif
dqRNGkind("Xoroshiro128+")
dqset.seed(1)
# RNG <- runif
# set.seed(1)

# prod.mean samples prodcnt U(0, 1) products of length N and
# calculates a summary.
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

# Uniform Product Distribution
# mathworld.wolfram.com/UniformProductDistribution.html
# u = x1 x ... xn
# density f(u, n) =abs(log(u))^(n-1) / (n-1)!
# integral (0, 1) f(u, n) * u du = 2^-n
# integral (0, 1) f(u, n) * log(u) du = n

# lambda = -log(value)
# dgamma(x=lambda, shape=k)   = dpois(k-1, lambda=lambda)
# 1-pgamma(x=lambda, shape=k) = ppois(k-1, lambda=lambda)


# WolframAlpha
# integrate x^(k-1) * e^-(2*x)  dx from 0 to inf = 2^(-k) Tau(k, 2)
# 2^-k * Tau(k, 0) / (k-1)! =  2^-k

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

prod.plotSampleMean <- function() {
    xmax <- 100
    xmin <- 50
    ymin <- 0
    ymax <- 100

    r1 <- c(1e+07, 50, 8.3e-16, 35, 50.1, 0.695, 1.002)
    r2 <- c(1e+07, 100, 1.4e-31, 71, 102.5, 0.711, 1.025)
    r3 <- c(1e+07, 150, 9.1e-48, 108, 156.3, 0.722, 1.042)
    r4 <- c(1e+07, 200, 8.2e-65, 148, 212.9, 0.738, 1.064)
    r5 <- c(1e+07, 250, 2.0e-83, 190, 274.7, 0.762, 1.099)
    r6 <- c(1e+07, 300, 1.8e-101, 232, 334.7, 0.773, 1.116)
    r7 <- c(1e+07, 400, 8.7e-141, 323, 465.3, 0.806, 1.163)
    r8 <- c(1e+07, 500, 8.0e-178, 408, 588.3, 0.816, 1.177)
    r9 <- c(1e+07, 750, 8.9e-273, 626, 903.7, 0.835, 1.205)
    r10 <- c(1e+07, 1000, 1.6e-367, 845, 1218.5, 0.845, 1.218)
    r11 <- c(1e+07, 1250, 1.1e-475, 1094, 1577.8, 0.875, 1.262)
    r12 <- c(1e+07, 1500, 3.8e-575, 1323, 1908.2, 0.882, 1.272)
    r13 <- c(1e+06, 2000, 3.5e-790, 1818, 2622.5, 0.909, 1.311)
    r14 <- c(1e+06, 2500, 9.7e-996, 2291, 3305.4, 0.916, 1.322)
    r15 <- c(1e+06, 3000, 6.0e-1195, 2750, 3967.1, 0.917, 1.322)
    r16 <- c(1e+06, 3500, 1.5e-1408, 3242, 4676.6, 0.926, 1.336)
    r17 <- c(1e+06, 4000, 2.5e-1614, 3715, 5360.3, 0.929, 1.340)
    r18 <- c(1e+06, 4500, 2.0e-1824, 4199, 6058.2, 0.933, 1.346)
    r19 <- c(1e+06, 5000, 1.5e-2040, 4697, 6776.2, 0.939, 1.355)
    r20 <- c(1e+06, 5500, 6.1e-2237, 5149, 7428.6, 0.936, 1.351)
    r21 <- c(1e+06, 6000, 1.2e-2455, 5653, 8155.1, 0.942, 1.359)
    r22 <- c(1e+06, 6500, 5.4e-2659, 6121, 8830.6, 0.942, 1.359)
    r23 <- c(1e+06, 7000, 1.4e-2879, 6629, 9563.3, 0.947, 1.366)
    r24 <- c(1e+06, 7500, 2.4e-3074, 7077, 10210.3, 0.944, 1.361)
    r25 <- c(1e+06, 8000, 1.4e-3302, 7603, 10968.6, 0.950, 1.371)
    r26 <- c(1e+06, 8500, 7.1e-3506, 8071, 11643.8, 0.950, 1.370)
    r27 <- c(1e+06, 9000, 1.8e-3728, 8583, 12383.3, 0.954, 1.376)
    r28 <- c(1e+06, 9500, 3.1e-3934, 9057, 13066.8, 0.953, 1.375)
    r29 <- c(1e+06, 10000, 1.9e-4132, 9514, 13725.3, 0.951, 1.373)
    r30 <- c(1e+06, 10500, 1.9e-4368, 10057, 14509.2, 0.958, 1.382)
    r31 <- c(1e+06, 11000, 2.3e-4570, 10522, 15180.0, 0.957, 1.380)

    mat <- rbind(
        r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13,
        r14, r15, r16, r17, r18, r19, r20, r21, r22, r23, r24,
        r25, r26, r27, r28, r29, r30, r31
    )
    cx <- mat[, 2]
    cy <- -100 * (mat[, 2] - mat[, 4]) / mat[, 2]
    cy2 <- 100 * (mat[, 2] - mat[, 5]) / mat[, 2]

    xlabels <- as.character(mat[, 2])
    plot.new()
    plot(cx, cy,
        xaxt = "n", # yaxt = "n",
        cex.main = 1, cex = 1, font.main = 1, ,
        col = "blue",
        ylim = c(-40, 0),

        xlab = "Number of variables in the product",
        main = "Sample mean difference (%) to e^N and 2^-N",
        ylab = "Difference (%) of logarithms"
    )
    par(new = TRUE)
    plot(cx, cy2,
        xaxt = "n", yaxt = "n",
        xlab = "", ylab = "",
        ylim = c(-40, 0), pch = 5,
        col = "red",
    )
    abline(h = 0, lty = 0.25, col = "#817e7e")
    abline(h = -4, lty = 2, lwd = 0.5, col = "blue")
    axis(1, labes <- xlabels, las = 2, tick = F)
    lines(cx, cy, lwd = 0.5, col = "blue")
    lines(cx, cy2, lwd = 0.5, col = "red")
}