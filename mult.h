#ifndef SpMV_H
#define SpMV_H

#include <cstdlib>
#include <cstdio>

typedef double number;


void SpMV(number *dst, number *src, number *result, number *vector, number *MAT, number *A, const int intervals);


#endif
