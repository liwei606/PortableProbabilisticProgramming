#include "ppp.h"

int flip_example_accept(int n, struct BNVertex* vertices[])  // This can be used as `accept`
{
    float x = vertices[5]->sample;  // vertices index starts from 0
    return (x > 2.0f && x < 3.0f);
}

int main()
{
    struct BNVertexDrawBern flip;
    struct BNVertexDrawNorm x1;
    struct BNVertexDrawGamma x2;
    struct BNVertexDrawConst $1;
    struct BNVertexComputePlus $2;
    struct BNVertexComputeIf x;

    struct BNVertex* verticies[6] = {
        (struct BNVertex*)&flip,
        (struct BNVertex*)&x1,
        (struct BNVertex*)&x2,
        (struct BNVertex*)&$1,
        (struct BNVertex*)&$2,
        (struct BNVertex*)&x
    };

    flip.super.super.type = DRAW;
    flip.super.type = FLIP;
    flip.p = 0.5;

    x1.super.super.type = DRAW;
    x1.super.type = GAUSSIAN;
    x1.mean = 0;
    x1.variance = 1;

    x2.super.super.type = DRAW;
    x2.super.type = GAMMA;
    x2.a = 1;
    x2.b = 1;

    $1.super.super.type = DRAW;
    $1.super.type = CONST;
    $1.c = 2;

    $2.super.super.type = COMPU;
    $2.super.type = PLUS;
    $2.left = (struct BNVertex*)&x2;
    $2.right = (struct BNVertex*)&$1;

    x.super.super.type = COMPU;
    x.super.type = IF;
    x.condition = (struct BNVertex*)&flip;
    x.consequent = (struct BNVertex*)&x1;
    x.alternative = (struct BNVertex*)&$2;

    rejection_sampling(6, verticies, flip_example_accept, print_vertices);
    return 0;
}
