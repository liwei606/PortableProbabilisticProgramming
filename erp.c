/*
Author: Wei Li

The xrps include:
bernoulli
log-flip 	            (A version of flip that takes log-probability parameter, to avoid underflow.)
sample-integer 	        (Sample a random integer between 0 and arg-1.)
sample-discrete         (Sample an integer according to a discrete distribution.)
gaussian
poisson
beta
gamma
exponential
dirichlet
uniform 	            (Uniform real number within a range.)
uniform-discrete        (Uniform integers within a range.)
no-proposals 	(this is the identity function, but tells the inference engine not to make proposals to any erp inside it. use this carefully -- it can save work when there is an erp that is known to be unchangeable, such as (equal? (erp) data), but can make the sampler wrong if used on something that should be sampled.)

*/

#include "stdio.h"
#include "stdlib.h"
#include "math.h"

// Generate 0 <= x <= 1;
double randomC() {
    return (double)rand() / RAND_MAX;
}

// Generate 0 < x <= 1;
double randomR() {
    return (long)(rand() + 1) / (long)(RAND_MAX + 1);
}

// Generate 0 <= x < 1;
double randomL() {
    return (long)rand() / (long)(RAND_MAX + 1);
}

int bernoulli(double p) {
	if (randomL() < p) return 1;
	else return 0;
}

int bernoulliD() {
    return bernoulli(0.5);
}

double uniform(double low, double high) {
    double v = randomC();
    return (1 - v) * low + v * high;
}

int uniformDiscrete(int low, int high) {
    double v = randomL();
    return (int)((1 - v) * low + v * (high + 1));
}

double gaussian(double mu, double sigma) {
    double u, v, x, y, q;
        do {
            u = 1 - randomC();
            v = 1.7156 * (randomC() - 0.5);
            x = u - 0.449871;
            y = fabs(v) + 0.386595;
            q = x * x + y * (0.196 * y - 0.25472 * x);
           }
        while (q >= 0.27597 && (q > 0.27846 || v * v > -4 * u * u * log(u)));
        return mu + sigma * v / u;
}

double gamma(double a, double b) {
    if (a < 1) return gamma(1 + a, b) * pow(randomC(), 1 / a);
        double x, v, u;
        double d = a - 1/3;
        double c = 1 / sqrt(9 * d);
        while (true)
            do {
                x = gaussian(0, 1);
                v = 1 + c * x;
               }
            while (v <= 0);
            v = v * v * v;
            u = randomC();
            if ((u < 1 - 0.331 * x * x * x * x) || (log(u) < 0.5 * x * x + d * (1 - v + log(v)))) 
                return b * d * v;
}

int main() {
    int i = 0;
    for (; i < 3000; ++i) {
        printf("%f\n", gaussian(0, 1));      
    }
    return 0;
}

