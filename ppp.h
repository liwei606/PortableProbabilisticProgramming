#ifndef PPP_H
#define PPP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

/************************************** Data Structures of the Distribution ******************************************/

enum {
    DRAW = 1,
    COMPU = 2
};

enum {
    IF = 1,
    PLUS = 2
};

enum {
    FLIP, LOG_FLIP,
    MULTINOMIAL,
    UNIFORM,
    UNIFORM_DISCRETE,
    GAUSSIAN,
    GAMMA,
    BETA,
    BINOMIAL,
    POISSON,
    DIRICHLET,
    CONST,
    DISTR_TYPE_LIMIT
} type;


struct BNVertex {
    int type;      // draw or compute
    float sample;  // last sampled value
};

struct BNVertexDraw {
    struct BNVertex super;  // extends BNVertex
    int type;               // dbern, dnorm, dgamma, or constant
};

struct BNVertexDrawBern {
    struct BNVertexDraw super;
    float p;
};

struct BNVertexDrawNorm {
    struct BNVertexDraw super;
    float mean;
    float variance;
};

struct BNVertexDrawGamma {
    struct BNVertexDraw super;
    float a;
    float b;
};

struct BNVertexDrawConst {
    struct BNVertexDraw super;
    float c;
};


struct BNVertexDrawBinomial {
	struct BNVertexDraw super;
	float p;
        int n;
};

struct BNVertexDrawMultinomial {
        struct BNVertexDraw super;
        float* theta;
        int n;
};

struct BNVertexDrawUniform {
        struct BNVertexDraw super;
        float low;
        float high;
};

struct BNVertexDrawUniformDiscrete {
        struct BNVertexDraw super;
        int low;
        int high;
};

struct BNVertexDrawBeta {
        struct BNVertexDraw super;
        float a;
        float b;
};

struct BNVertexDrawPoisson {
        struct BNVertexDraw super;
        int mu;
};

struct BNVertexDrawDrichlet {
        struct BNVertexDraw super;
        float* alpha;
        int n;
};

struct BNVertexCompute {
    struct BNVertex super;  // extends BNVertex
    int type;               // + or if
};


struct BNVertexComputePlus {
    struct BNVertexCompute super;
    struct BNVertex* left;
    struct BNVertex* right;
};

struct BNVertexComputeSubstract {
        struct BNVertexCompute super;
        struct BNVertex* left;
        struct BNVertex* right;
};

struct BNVertexComputeMultiply {
        struct BNVertexCompute super;
        struct BNVertex* left;
        struct BNVertex* right;
};

struct BNVertexComputeDivide {
        struct BNVertexCompute super;
        struct BNVertex* left;
        struct BNVertex* right;
};

struct BNVertexComputeIf {
    struct BNVertexCompute super;
    struct BNVertex* condition;
    struct BNVertex* consequent;
    struct BNVertex* alternative;
};

/************************************** Sample Functions ******************************************/

float flip(float p);
float flipD();
float log_flip(float p);

float multinomial(float theta[], int n);
float multinomial_logprob(float m, float theta[], int n);

float uniform(float low, float high);
float uniformDiscrete(int low, int high);

float gaussian(float mu, float sigma);
float guassian_logprob(float x, float mu, float sigma);

float gamma(float a, float b);
float gamma_logprob(float x, float a, float b);

float beta(float a, float b);
float beta_logprob(float x, float a, float b);

float binomial(float p, int n);
float binomial_lgoprob(float s, float p, int n);

float poisson(int mu);
float poisson_logprob(float k, int mu);

float* dirichlet(float alpha[], int n);
float dirchlet_lgoprob(float theta[], float alpha[]);

/************************************** Inference Functions ******************************************/

typedef int (*vertices_handler)(int n, struct BNVertex* vertices[]);
typedef int (*inference_engine)(int n, struct BNVertex* vertices[], vertices_handler accept, vertices_handler add_trace);

int rejection_sampling(int n, struct BNVertex* vertices[], vertices_handler accept, vertices_handler add_trace);

int print_vertices(int n, struct BNVertex* vertices[]);

#endif
