# uniprod

Here are functions for computing and comparing
sample mean values of products of Uniform(0, 1)
random numbers.

Loading prod.R in workspace calls Rcpp/C++ compiler for C++ functions in uniprod.cpp file. Copy both files to a same directory.

R functions implement most of the functionality and C++ functions are not entirely necessary.

The functions are some what commented in the program code.

## Examples

If parameter gamma = TRUE, product of N U(0, 1) variables
is simulated by a single gamma random variable:

```R
    prod(runif(N)) <is distributed> exp(-rgamma(1,N))
```

```R
# Multiplies U(0, 1) variables until product is zero.
# Returns number of variables.
prod.tozero()

# R function, which uses R mean(vector) for sample mean.
# Probably gives the best estimate of product mean.
prod.mean(samplesize = 1e6, N = 200, gamma = F, seed = 0)

# C++ function with 80-bit extended precision computation
prod.mean80(samplesize = 1e7, N = 5000, gamma = T, seed = 0)

# Faster C++ version of prod.mean. Uses own summing control
# for the mean computation.
prod.mean64(samplesize = 1e7, N = 200, gamma = F, seed = 0)

# prod.long(N) multiplies N U(0, 1) random variables and
# compares the log(result) to log(e^-N) = N.
prod.long(N = 1e7)

# Runs samplesize products to e^lim and plots sampled product
# lenght distribution against Poisson(-lim) distribution.
prod.plotPDF(lim = -745, samplesize = 1e+6, log = T)
```
