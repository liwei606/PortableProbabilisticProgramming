#ifndef PPP_H
#define PPP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

/************************************** Data Structures of the Distribution ******************************************/

#define DRAW 1
#define COMPU 2

#define IF 1
#define PLUS 2

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
 
struct BNVertexCompute {
  struct BNVertex super;  // extends BNVertex
  int type;               // + or if
};
 
struct BNVertexComputePlus {
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

typedef int (*vertices_handler)(int n, struct BNVertex* vertices[]);
typedef int (*inference_engine)(int n, struct BNVertex* vertices[], vertices_handler accept, vertices_handler add_trace);


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


#endif
