
#include <Rcpp.h>
using namespace Rcpp;

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}
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
// [[Rcpp::export]]
int setseed(uint64_t seed) {
    uint64_t z = seed + 0x9e3779b97f4a7c15;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    s0 = z ^ (z >> 31);
    s1 = z + 0x9e3779b97f4a7c15;
    return 0;
}

// prodrunif is vectorized prng based on xoroshiro128+.
// [[Rcpp::export]]
NumericVector prodrunif(int N) {
    NumericVector r(N);
    for (int i = 0; i < N; i++) {
        r[i] = myrunif();
    }
    return r;
}

// This returns exponent of x: 0 - 1022 for x < 1.
static inline int exponent(double x) {
    const union {
        uint64_t u;
        double d;
    } z = {.d = x};
    return (int)(z.u >> 52) & 0x7ff;
}
// prodmean64 uses array mexp[1023] to sum products of same magnitude
// to a same array slot.
// [[Rcpp::export]]
int prodmean64(double samplesize = 1e+6, double N = 100,
               bool gamma = false, uint64_t seed = 0) {
    double meanprod = 0, meanlog = 0, logp = 0, p, pmax = 0;
    double mexp[1023] = {0};
    if (!gamma && seed > 0) setseed(seed);

    for (uint64_t i = 0; i < samplesize; i++) {
        if (gamma) {
            logp = -R::rgamma(N, 1);
            p = exp(logp);
        } else {
            p = 1;
            for (int k = 0; k < N; k++)
                p *= myrunif();
            if (p > 0) logp = log(p);
        }
        mexp[exponent(p)] += p;
        meanlog += logp;
        if (p > pmax) {
            pmax = p;
        }
    }
    for (int i = 0; i < 1023; i++) {
        meanprod += mexp[i]; // from small to big
    }
    meanprod /= samplesize;
    meanlog /= samplesize;

    printf("sample size       %1.0e\n", samplesize);
    printf("N                 %1.0f\n", N);
    printf("2^-N              %1.0e  (mean distribution)\n", pow(2, -N));
    printf("mean sample       %1.0e  %1.15e\n", meanprod, meanprod);
    printf("e^-N              %1.0e\n", exp(-N));
    printf("max sample        %1.0e\n", pmax);
    printf("log mean          %1.0f\n", log(meanprod));
    printf("mean log(prod)    %1.0f\n", meanlog);
    return 0;
}

// Extended precision by long doubles gives 19 decimal digits accurary
// and 16 x larger exponent space, which is more relevant here.
// prod.mean80(1e6, N=100, gamma=F, seed=0)
// [[Rcpp::export]]
int prodmean80(long double samplesize = 1e+6, long double N = 100,
               bool gamma = false, uint64_t seed = 0) {
    long double meanprod = 0, meanlog = 0,
                logp = 0, nextsum = 0, sum = 0, missed = 0, p;
    if (!gamma && seed > 0) setseed(seed);

    for (uint64_t i = 0; i < samplesize; i++) {
        if (gamma) {
            logp = -R::rgamma(N, 1);
            p = expl(logp);
        } else {
            p = 1;
            for (int k = 0; k < N; k++)
                p *= myrunif();
            if (p > 0) logp = log(p);
        }
        meanlog += logp;
        // Kahan/Neumaier summation
        nextsum = sum + p;
        if (sum >= p)
            missed += (sum - nextsum) + p;
        else
            missed += (p - nextsum) + sum;
        sum = nextsum;
    }
    sum += missed;
    meanprod = sum / samplesize;
    meanlog /= samplesize;

    printf("sample size       %1.0Le\n", samplesize);
    printf("N                 %1.0Lf\n", N);
    printf("2^-N              %1.0Le\n", powl(2, -N));
    printf("mean sample       %1.0Le  %1.18Le\n", meanprod, meanprod);
    printf("e^-N              %1.0Le\n", expl(-N));
    printf("log mean          %1.0Lf\n", logl(meanprod));
    printf("mean log(prod)    %1.0Lf\n", meanlog);
    return 0;
}

#include <quadmath.h>

// prodmean128 uses 128-bit quadruple-precision arithmetic for mean.
// Quadruple-precision gives 34 decimal digits accuracy.
// [[Rcpp::export]]
int prodmean128(long double samplesize = 1e+6, long double N = 100, uint64_t seed = 0) {
    long double meanprod;
    __float128 sum = 0.0, p;
    if (seed > 0) setseed(seed);

    for (uint64_t i = 0; i < samplesize; i++) {
        p = 1.0;
        for (int k = 0; k < N; k++)
            p *= (__float128)myrunif();
        sum += p;
    }
    meanprod = (long double)(sum / (__float128)samplesize);

    printf("sample size       %1.0Le\n", samplesize);
    printf("N                 %1.0Lf\n", N);
    printf("2^-N              %1.0Le\n", powl(2, -N));
    printf("mean sample       %1.0Le  %1.18Le\n", meanprod, meanprod);
    printf("e^-N              %1.0Le\n", expl(-N));
    return 0;
}
