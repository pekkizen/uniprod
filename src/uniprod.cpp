
#include <cmath>
#include <cstdio>
#include <stdint.h>

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

// prng.di.unimi.it/xoroshiro128plus.c
static uint64_t s[2] = {0x83b5b142866da9d5, 0xcb9c59b3f9f87d4d};
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

void samplemean(int rounds, int N, int call) {
    long double prod, lgm, lg2m, mean = 0;

    for (int i = 0; i < rounds; i++) {
        prod = 1;
        for (int k = 0; k < N; k++) {
            prod *= rand64();
        }
        mean += prod;
    }
    mean /= rounds;
    lgm = logl(mean);
    lg2m = log2l(mean);
    // printf("r%d <- c(%1.0e, %d, %1.1Le, %1.0Lf, %1.1Lf, %1.3Lf, %1.3Lf)\n",
    //        call, (double)rounds, N, mean, -lgm, -lg2m, -lgm / N, -lg2m / N);

    printf("rounds        %d\n", rounds);
    printf("N             %d\n", N);
    printf("mean          %1.1Le\n", mean);
    printf("log(mean)     %1.0Lf\n", logl(mean));
    printf("log2(mean)    %1.0Lf\n", log2l(mean));
    printf("log(mean)/ N  %1.4Lf\n", -logl(mean) / N);
}

// uniprod maxN, rounds, seed
// uniprod 10000, 1e4, 123
int main(int an, char **arg) {

    if (an < 4) return 1;

    long N = atol(arg[1]);
    int rounds = atol(arg[2]);
    long seed = atol(arg[3]);
    if (seed > 0) {
        s[0] *= seed;
        s[1] += seed;
    }
    if (N > 11000) return 1;
    int call = 0;
    int increment = 50;
    for (int k = 50; k <= N; k += increment) {
        call++;
        samplemean(rounds, k, call);
        if (k >= 300) increment = 100;
        if (k >= 500) increment = 250;
        if (k >= 2000) increment = 500;
    }
    return 0;
}
