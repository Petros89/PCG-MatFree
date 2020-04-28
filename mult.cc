#include <cstdlib>
#include <cstdio>
#include "mult.h"

typedef double number;

void SpMV(number *dst, number *src, number *result, number *vector, number *MAT, number *A, const int intervals){

  const int dim = 3;
  const int intervals_in_color = intervals >> 1;
  const int N = intervals + 1;
  const int ndofs = 1<<dim;

  // loop over all colors
  for (int color3=0; color3 <2; color3++) {
    for (int color2=0; color2 <2; color2++) {
      for (int color1=0; color1 <2; color1++) {


        // loop over all elements
        for (int k3=0; k3 < intervals_in_color; k3++) {
          for (int k2=0; k2 < intervals_in_color; k2++) {
            for (int k1=0; k1 < intervals_in_color; k1++) {

              const int K = (color1 + k1*2) + intervals*((color2 + k2*2) + intervals*(color3 + k3*2));

              //---------------------------------------------------------------
              // element-local operation
              //---------------------------------------------------------------

              // read local dof values
              for (int i3=0; i3 <2; i3++) {
                for (int i2=0; i2 <2; i2++) {
                  for (int i1=0; i1 <2; i1++) {
                    const int global_idx = (i1 + color1 + k1*2) + N*((i2 + color2 + k2*2) + N*(i3 + color3 + k3*2));
                    result[i1+2*(i2+2*i3)] = src[global_idx];
                  }
                }
              }


              number tmp = 0;
              for(int i=0; i<ndofs; i++) {

                for(int j=0; j<ndofs; j++) {
                   tmp += MAT[i*ndofs+j]*result[j];
                }
		A[K]=1.0;
    	        vector[i] = tmp*A[K];
              }



              // write back local dof values
              for (int i3=0; i3 <2; i3++) {
                for (int i2=0; i2 <2; i2++) {
                  for (int i1=0; i1 <2; i1++) {
                    const int global_idx = (i1 + color1 + k1*2) + N*((i2 + color2 + k2*2) + N*(i3 + color3 + k3*2));
                    dst[global_idx] += vector[i1+2*(i2+2*i3)];
                  }
                }
              }

            }
          }
        }


      }
    }
  }

}




