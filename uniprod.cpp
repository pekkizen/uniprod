
#include <cmath>
#include <cstdio>
#include <stdint.h>

#include <Rcpp.h>
using namespace Rcpp;

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}
// Random state for xoroshiro128plus
static uint64_t s0 = 0x9e3779b97f4a7c15;
static uint64_t s1 = 0xbf58476d1ce4e5b9;

// prng.di.unimi.it
static inline uint64_t xoroshiro128plus(void) {
    const uint64_t result = s0 + s1;
    s1 ^= s0;
    s0 = rotl(s0, 24) ^ s1 ^ (s1 << 16);
    s1 = rotl(s1, 37);
    return result;
}
static inline double myrunif() {
    return ((double)(xoroshiro128plus() >> 11)) * 0x1.0p-53;
}

// Adapted from prng.di.unimi.it/splitmix64.c
// [[Rcpp::export(name = prod.set.seed)]]
int setseed(uint64_t seed) {
    uint64_t z = seed + 0x9e3779b97f4a7c15;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    s0 = z ^ (z >> 31);
    s1 = z + 0x9e3779b97f4a7c15;
    return 0;
}

// [[Rcpp::export(name = prod.runif)]]
NumericVector prodrunif(int N) {
    NumericVector r(N);
    for (int i = 0; i < N; i++) {
        r[i] = myrunif();
    }
    return r;
}

// This gets exponent of x: 0 - 1022 for x < 1.
static inline int exponent(double x) {
    const union {
        uint64_t u;
        double d;
    } z = {.d = x};
    return (int)(z.u >> 52) & 0x7ff;
}

// https://en.wikipedia.org/wiki/Extended_precision
// Extended precision gives 3 more digits accurary
// and 16 x larger exponent space, which is more relevant here.
// prod.mean80(1e6, N=100, gamma=F, seed=0)
// [[Rcpp::export]]
int prodmean80(long double samplesize = 1e+6, long double N = 100,
               bool gamma = false, uint64_t seed = 0) {
    long double meanprod = 0, meanlog = 0, meangeom = 0, logprod = 0,
                prev = 0, missedsum = 0, prod;
    if (!gamma && seed > 0) setseed(seed);

    for (long i = 0; i < samplesize; i++) {
        if (gamma) {
            logprod = R::rgamma(N, 1);
            prod = expl(-logprod);
        } else {
            prod = 1;
            for (int k = 0; k < N; k++)
                prod *= myrunif();
        }
        meanprod += prod;
        if (meanprod == prev) missedsum += prod;
        prev = meanprod;
        meanlog += logprod;
        meangeom += expl(logl(prod) / N);
    }
    meanprod += missedsum; // doesn't help
    meanprod /= samplesize;
    meanlog /= samplesize;
    meangeom /= samplesize;

    printf("sample size       %1.0Le\n", samplesize);
    printf("2^-N              %1.0Le\n", powl(2, -N));
    printf("mean(prod)        %1.0Le\n", meanprod);
    printf("e^-N              %1.0Le\n", expl(-N));
    printf("log(mean)         %1.0Lf\n", logl(meanprod));
    printf("mean(log(prod))   %1.0Lf\n", -meanlog);
    printf("mean geom         %1.5Lf\n", meangeom);
    printf("(1 + 1/N)^-N      %1.5Lf\n", powl(1 + 1 / N, -N));
    printf("1/e               %1.5f\n", 1 / exp(1));
    printf("2^-N/mean         %1.0Le\n", powl(2, -N) / meanprod);
    printf("e^-N/mean         %1.0Le\n", expl(-N) / meanprod);
    printf("mean              %1.18Le\n", meanprod);
    printf("mean missed sum   %1.0Le\n", missedsum / samplesize);
    return 0;
}

// prodmean64 uses array mexp[1023] to sum product of same magnitude
// to a same array slot. mexp[1023] represents all possible exponent
// values for 64-bit floats < 1. Number x is added to slot[e], where
// e = exponent of x: mexp[exponent(prod)] += prod;
//
// [[Rcpp::export]]
int prodmean64(double samplesize = 1e+6, double N = 100,
               bool gamma = false, uint64_t seed = 0) {
    double meanprod = 0, meanlog = 0, meangeom = 0,
           logprod = 0, prod;
    double mexp[1023] = {0};
    if (!gamma && seed > 0) setseed(seed);

    for (long i = 0; i < samplesize; i++) {
        if (gamma) {
            logprod = R::rgamma(N, 1);
            prod = exp(-logprod);
        } else {
            prod = 1;
            for (int k = 0; k < N; k++)
                prod *= myrunif();
            if (prod > 0) logprod = -logl(prod);
        }
        mexp[exponent(prod)] += prod;

        meanlog += logprod;
        meangeom += exp(log(prod) / N);
    }
    for (int i = 0; i < 1023; i++) {
        meanprod += mexp[i]; // from small to big
    }
    meanprod /= samplesize;
    meanlog /= samplesize;
    meangeom /= samplesize;

    printf("sample size       %1.0e\n", samplesize);
    printf("2^-N              %1.0e\n", pow(2, -N));
    printf("mean(prod)        %1.0e\n", meanprod);
    printf("e^-N              %1.0e\n", exp(-N));
    printf("log(mean)         %1.0f\n", log(meanprod));
    printf("mean(log(prod))   %1.0f\n", -meanlog);
    printf("mean geom         %1.5f\n", meangeom);
    printf("(1 + 1/N)^-N      %1.5f\n", pow(1 + 1 / N, -N));
    printf("1/e               %1.5f\n", 1 / exp(1));

    printf("2^-N/mean         %1.0e\n", pow(2, -N) / meanprod);
    printf("e^-N/mean         %1.0e\n", exp(-N) / meanprod);
    printf("mean              %1.16e\n", meanprod);
    return 0;
}
