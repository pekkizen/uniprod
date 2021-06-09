# uniprod

Here are functions for computing and comparing
sample mean values of products of Uniform(0, 1)
random numbers.

Loading prod.R in workspace calls Rcpp/C++ compiler for C++ functions in uniprod.cpp file. Copy both files to a same directory.

R functions implement most of the functionality and C++ functions are not entirely necessary.

The functions are some what commented in the program code.

## Examples

```R
prod.tozero()

prod.mean(samplesize = 1e6, N = 200, gamma = F, seed = 0)
prod.mean80(samplesize = 1e7, N = 5000, gamma = T, seed = 0)
prod.mean64(samplesize = 1e7, N = 200, gamma = F, seed = 0)

prod.long(N = 1e7)
prod.plotPDF(lim = -745, samplesize = 1e+6, log = T)
```
