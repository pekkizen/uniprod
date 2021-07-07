# uniprod

Functions for computing statistics of products of Uniform(0, 1) random numbers.

Loading prod.R in workspace calls Rcpp/C++ compiler for C++ functions in uniprod.cpp file. Copy both files to a same directory.

R functions implement most of the functionality and C++ functions are not entirely necessary.

## Examples

If parameter gamma = TRUE, the product of N U(0, 1) variables
is simulated by a single Gamma random variable:

```R
    prod(runif(N)) is simulated by exp(-rgamma(1,N))
```

If parameter seed > 0, random number generator is seeded by the seed.

```R
# Multiplies U(0, 1) variables until product is zero.
# Returns number of variables.
prod.tozero()

# prod.mean calculates sample mean and median and statistics.
prod.mean(samplesize = 1e6, N = 200, gamma = F, seed = 0)

# Function for calculating sample mean
# C++ function with 80-bit extended precision computation
# Can handle products up to 11000 variables.
prod.mean80(samplesize = 1e7, N = 5000, gamma = T, seed = 0)

# prod.clim calculates and tests confidence intervals.
prod.clim(samplesize = 1e6, N = 200, gamma = F, seed = 0)

# prod.long(N) multiplies N U(0, 1) random variables and
# compares the log(result) to log(e^-N) = N.
prod.long(N = 1e7)

# Runs samplesize products to e^lim and plots sampled product
# lenght distribution against Poisson(-lim) distribution.
prod.plotPDF(lim = -500, samplesize = 1e+5, log = T)
```
