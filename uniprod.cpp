
#include <cmath>
#include <cstdio>
#include <stdint.h>

#include <Rcpp.h>
using namespace Rcpp;

// https://en.wikipedia.org/wiki/Extended_precision

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}
static uint64_t s[2] = {0x9e3779b97f4a7c15, 0xbf58476d1ce4e5b9};

// prng.di.unimi.it/xoroshiro128plus.c
static inline uint64_t xoroshiro128plus(void) {
    const uint64_t s0 = s[0];
    uint64_t s1 = s[1];
    const uint64_t result = s0 + s1;
    s1 ^= s0;
    s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16);
    s[1] = rotl(s1, 37);
    return result;
}
static inline double rand64() {
    return ((double)(xoroshiro128plus() >> 11)) * 0x1.0p-53;
}

// [[Rcpp::export(name = prod.set.seed)]]
int setseed(uint64_t seed) {
    uint64_t z = seed + 0x9e3779b97f4a7c15;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    s[0] = z ^ (z >> 31);
    s[1] = z + 0xf7c2ebc08f67f2b5;
    return 0;
}

// [[Rcpp::export(name = prod.runif)]]
NumericVector prodrunif(int N) {
    NumericVector r(N);
    for (int i = 0; i < N; i++) {
        r[i] = rand64();
    }
    return r;
}

// [[Rcpp::export(name = prod.mean80)]]
int prodmean80(long double samplesize = 1e+6, long double N = 100) {
    long double meanprod = 0, meanlog = 0,
                meangeom = 0, logprod = 0, prod;

    for (long i = 0; i < samplesize; i++) {
        if (N <= 300) {
            prod = 1;
            for (int k = 0; k < N; k++)
                prod *= rand64();
            logprod = -logl(prod);
        } else {
            logprod = (long double)R::rgamma(N, 1);
            prod = expl(-logprod);
        }
        meanprod += prod;
        meanlog += logprod;
        meangeom += powl(prod, 1.0 / N);
    }
    meanprod /= samplesize;
    meanlog /= samplesize;
    meangeom /= samplesize;

    printf("sample size       %1.0Le\n", samplesize);
    printf("mean              %1.1Le\n", meanprod);
    printf("log(mean)         %1.1Lf\n", -logl(meanprod));
    printf("mean(log(prod))   %1.1Lf\n", meanlog);
    printf("mean geom         %1.5Lf\n", meangeom);
    printf("(1 + 1/N)^-N      %1.5Lf\n", powl(1 + 1 / N, -N));
    printf("1/e               %1.5f\n", 1 / exp(1));
    printf("2^-N              %1.0Le\n", powl(2, -N));
    printf("2^-N/mean         %1.0Le\n", powl(2, -N) / meanprod);
    printf("e^-N/mean         %1.0Le\n", expl(-N) / meanprod);
    printf("mean              %1.19Le\n", meanprod);
    return 0;
}

static inline int exponent(double x) {
    const union {
        uint64_t i;
        double d;
    } u = {.d = x};
    return (int)(u.i >> 52) & 0x7ff;
}

// [[Rcpp::export(name = prod.mean64)]]
int prodmean64(double samplesize = 1e+6, double N = 100) {
    double meanprod = 0, meanlog = 0,
           meangeom = 0, logprod = 0, prod;

    double mexp[1023] = {0};
    for (long i = 0; i < samplesize; i++) {

        if (samplesize <= 1e7) {
            prod = 1;
            for (int k = 0; k < N; k++)
                prod *= rand64();
            logprod = -log(prod);
        } else {
            logprod = R::rgamma(N, 1);
            prod = exp(-logprod);
        }
        mexp[exponent(prod)] += prod;

        meanlog += logprod;
        meangeom += pow(prod, 1.0 / N);
    }
    for (int i = 0; i < 1023; i++)
        meanprod += mexp[i]; // from small to big

    meanprod /= samplesize;
    meanlog /= samplesize;
    meangeom /= samplesize;

    printf("sample size       %1.0e\n", samplesize);
    printf("mean              %1.1e\n", meanprod);
    printf("log(mean)         %1.1f\n", log(meanprod));
    printf("mean(log(prod))   %1.1f\n", -meanlog);
    printf("mean geom         %1.5f\n", meangeom);
    printf("(1 + 1/N)^-N      %1.5f\n", pow(1 + 1 / N, -N));
    printf("1/e               %1.5f\n", 1 / exp(1));
    printf("2^-N              %1.0e\n", pow(2, -N));
    printf("2^-N/mean         %1.0e\n", pow(2, -N) / meanprod);
    printf("e^-N/mean         %1.0e\n", exp(-N) / meanprod);
    printf("mean              %1.16e\n", meanprod);
    return 0;
}
