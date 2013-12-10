/*

   Inference algorithms:
   1. Sampling
   Rejection Sampling
   Metropolis-Hastings Algorithm
   2. Exact Inference

 */

#include "ppp.h"

float getsample(struct BNVertexDraw* vertexDraw) {
    int distribution = vertexDraw->type;
    switch (distribution) {
        case FLIP: return flip(((struct BNVertexDrawBern*)vertexDraw)->p);
        case LOG_FLIP: return log_flip(((struct BNVertexDrawBern*)vertexDraw)->p);
        case GAUSSIAN: return gaussian(((struct BNVertexDrawNorm*)vertexDraw)->mean, ((struct BNVertexDrawNorm*)vertexDraw)->variance);
        case GAMMA: return gamma(((struct BNVertexDrawGamma*)vertexDraw)->a, ((struct BNVertexDrawGamma*)vertexDraw)->b);
        case CONST: return ((struct BNVertexDrawConst*)vertexDraw)->c;
        default: return -1;
    }
}

float getcomp(struct BNVertexCompute* vertexComp) {
    int operator = vertexComp->type;
    switch (operator) {
        case PLUS:
            return ((struct BNVertexComputePlus*)vertexComp)->left->sample + ((struct BNVertexComputePlus*)vertexComp)->right->sample;
        case IF: {
                     if (((struct BNVertexComputeIf*)vertexComp)->condition->sample == 1.0)
                         return ((struct BNVertexComputeIf*)vertexComp)->consequent->sample;
                     return ((struct BNVertexComputeIf*)vertexComp)->alternative->sample;
                 }
        default: return -1;
    }
}

int rejection_sampling(int n, struct BNVertex* vertices[], vertices_handler accept, vertices_handler add_trace) {
    int k, i;
    for (k = 0; k < 1000; ++k) {
        for (i = 0; i < n; ++i) {
            if (vertices[i]->type == DRAW) {
                vertices[i]->sample = getsample((struct BNVertexDraw*)vertices[i]);
            }
            else {
                vertices[i]->sample = getcomp((struct BNVertexCompute*)vertices[i]);
            }
        }
        if (accept(n, vertices)) {
            add_trace(n, vertices);
        }
    }
    return 0;
}


int print_vertices(int n, struct BNVertex* vertices[])  // This can be used as `add_trace`
{
    int i;
    for (i = 0; i < n; i++) {
        printf("%f ", vertices[i]->sample);
    }
    printf("\n");
    return 0;
}

