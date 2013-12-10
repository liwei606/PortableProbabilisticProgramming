/*

   The xrps include:
   flip
   log-flip	            (A version of flip that takes log-probability parameter, to avoid underflow.)
   multinomial
   uniform 	            (Uniform real number within a range.)
   uniform-discrete        (Uniform integers within a range.)
   gaussian
   gamma
   beta
   binomial
   poisson
   dirichlet
   no-proposals 	(this is the identity function, but tells the inference engine not to make proposals to any erp inside it. use this carefully -- it can save work when there is an erp that is known to be unchangeable, such as (equal? (erp) data), but can make the sampler wrong if used on something that should be sampled.)

 */

#include "ppp.h"
#include <math.h>

/************************************** Base Random Functions ******************************************/

// Generate 0 <= x <= 1;
float randomC() {
    return (float)rand() / RAND_MAX;
}

// Generate 0 < x <= 1;
float randomR() {
    return (long double)(rand() + 1) / ((long double)RAND_MAX + 1);
}

// Generate 0 <= x < 1;
float randomL() {
    return (long double)rand() / ((long double)RAND_MAX + 1);
}

int Round(float number) {
    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
}

/************************************** Flip ******************************************/

float flip(float p) {
    if (randomL() < p) return 1;
    return 0;
}

float flipD() {
    return flip(0.5);
}

float log_flip(float p) {
    if (log(randomL()) < p) return 1;
    return 0;
}

/************************************** Multinomial ******************************************/

float multinomial(float theta[], int n) {
    int k = n;
    float thetasum = 0;
    for (int i = 0; i < k; i++) thetasum += theta[i];
    float x = randomC() * thetasum;
    float probAccum = 0;
    for(int i = 0; i < k; i++) {
        probAccum += theta[i];
        if(probAccum >= x) return i; //FIXME: if x = 0 returns i=0, but this isn't right if theta[0] == 0...
    }
    return k;
}

float multinomial_logprob(float m, float theta[], int n) {
    int k = n;
    if (m < 0 || m >= k) {
        return INT_MIN;
    }
    m = Round(m);
    float thetasum = 0;
    for (int i = 0; i < k; i++) thetasum += theta[i];
    return log(theta[(int)m] / thetasum);
}

/************************************** Uniform ******************************************/

float uniform(float low, float high) {
    float v = randomC();
    return (1 - v) * low + v * high;
}

float uniformDiscrete(int low, int high) {
    float v = randomL();
    return (int)((1 - v) * low + v * (high + 1));
}

/************************************** Gaussian ******************************************/

float gaussian(float mu, float sigma) {
    float u, v, x, y, q;
    do {
        u = 1 - randomC();
        v = 1.7156 * (randomC() - 0.5);
        x = u - 0.449871;
        y = fabs(v) + 0.386595;
        q = x * x + y * (0.196 * y - 0.25472 * x);
    } while (q >= 0.27597 && (q > 0.27846 || v * v > -4 * u * u * log(u)));
    return mu + sigma * v / u;
}

float gaussian_logprob(float x, float mu, float sigma) {
    return -0.5 * (1.8378770664093453 + 2 * log(sigma) + (x - mu) * (x - mu) / (sigma * sigma));
}

/************************************** Gamma ******************************************/

float gamma(float a, float b) {
    if (a < 1) return gamma(1 + a, b) * pow(randomC(), 1 / a);
    float x, v, u;
    float d = a - 1/3;
    float c = 1 / sqrt(9 * d);
    while (1) {
        do {
            x = gaussian(0, 1);
            v = 1 + c * x;
        } while (v <= 0);
        v = v * v * v;
        u = randomC();
        if ((u < 1 - 0.331 * x * x * x * x) || (log(u) < 0.5 * x * x + d * (1 - v + log(v))))
            return b * d * v;
    }
}

float gamma_cof[6] = {76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};

float log_gamma(float xx) {
    float x = xx - 1;
    float tmp = x + 5.5;
    tmp = tmp - (x + 0.5) * log(tmp);
    float ser = 1.000000000190015;
    for (int j = 0; j <= 5; ++j) {
        x = x + 1;
        ser = ser + gamma_cof[j] / x;
    }
    return -tmp + log(2.5066282746310005 * ser);
}

float gamma_logprob(float x, float a,float b) {
    return (a - 1) * log(x) - x / b - log_gamma(a) - a * log(b);
}

/************************************** Beta ******************************************/

float beta(float a, float b) {
    float x = gamma(a, 1);
    return x / (x + gamma(b, 1));
}

float log_beta(float a, float b) {
    return log_gamma(a) + log_gamma(b) - log_gamma(a+b);
}

float beta_logprob(float x, float a, float b) {
    if (x > 0 && x < 1)
        return (a - 1) * log(x) + (b - 1) * log(1 - x) - log_beta(a, b);
    return INT_MIN;
}

/************************************** Binomial ******************************************/

float binomial(float p, int n) {
    int k = 0;
    int N = 10;

    float a, b;
    while (n > N) {
        a = 1 + n / 2;
        b = 1 + n - a;

        float x = beta(a,b);
        if (x >= p) {
            n = a - 1;
            p /= x;
        }
        else {
            k += a;
            n = b - 1;
            p = (p - x) / (1 - x);
        }
    }
    float u;
    for (int i = 0; i < n; i++) {
        u = randomC();
        if (u < p) k++;
    }
    return k;
}

float g(float x) {
    if (x == 0) return 1;
    if (x == 1) return 0;
    float d = 1 - x;
    return (1 - (x * x) + (2 * x * log(x))) / (d * d);
}

float binomial_logprob(float s, float p, int n) {
    float inv2 = 1/2;
    float inv3 = 1/3;
    float inv6 = 1/6;
    if (s >= n) return INT_MIN;
    float q = 1 - p;
    float S = s + inv2;
    float T = n - s - inv2;
    float d1 = s + inv6 - (n + inv3) * p;
    float d2 = q / (s + inv2) - p / (T + inv2) + (q - inv2) / (n + 1);
    d2 = d1 + 0.02 * d2;
    float num = 1 + q * g(S / (n * p)) + p * g(T / (n * q));
    float den = (n + inv6) * p * q;
    float z = num / den;
    float invsd = sqrt(z);
    z = d2 * invsd;
    return gaussian_logprob(z, 0, 1) + log(invsd);
}

/************************************** Poisson ******************************************/

float poisson(int mu) {
    int k = 0;
    while (mu > 10) {
        int m = 7 / 8 * mu;
        float x = gamma(m, 1);
        if (x > mu) return k + binomial(mu / x, m - 1);
        else {
            mu -= x;
            k += m;
        }
    }
    float emu = exp(-mu);
    float p = 1;
    do {
        p *= randomC();
        k++;
    } while (p > emu);
    return k-1;
}

int fact(int x) {
    int t = 1;
    while (x > 1) t *= x--;
    return t;
}


int lnfact(float x) {
    if (x < 1) x = 1;
    if (x < 12) return log(fact(Round(x)));
    float invx = 1 / x;
    float invx2 = invx * invx;
    float invx3 = invx2 * invx;
    float invx5 = invx3 * invx2;
    float invx7 = invx5 * invx2;

    float sum = ((x + 0.5) * log(x)) - x;
    sum += log(2 * 3.141592654) / 2;
    sum += (invx / 12) - (invx3 / 360);
    sum += (invx5 / 1260) - (invx7 / 1680);
    return sum;
}

float poisson_logprob(float k, int mu) {
    return k * log(mu) - mu - lnfact(k);
}

/************************************** Dirichlet ******************************************/

float* dirichlet(float alpha[], int n) {
    float ssum = 0;
    float *theta;
    theta = malloc(sizeof(float) * n);
    float t;
    for (int i = 0; i < n; i++) {
        t = gamma(alpha[i], 1);
        theta[i] = t;
        ssum = ssum + t;
    }
    for (int i = 0; i < n; i++) theta[i] /= ssum;
    return theta;
}

float dirichlet_logprob(float theta[], float alpha[], int n) {
    float asum = 0;
    for (int i = 0; i < n; i++) asum += alpha[i];
    float logp = log_gamma(asum);
    for (int i = 0; i < n; i++) {
        logp += (alpha[i] - 1) * log(theta[i]);
        logp -= log_gamma(alpha[i]);
    }
    return logp;
}

void test() {
    srand(time(NULL));
    printf("%f\n", randomC());
    printf("%f\n", randomL());
    printf("%f\n", flip(0.5));
}

