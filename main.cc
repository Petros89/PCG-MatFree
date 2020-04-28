/* -*- Matrix-Free Preconditioned Conjugate Gradient (PCG) -*-
 *
 * @(#)main.cc --- 01/27/2020
 * @author Petos Aposotolou <pea11@pitt.edu>
 *
 */

#include <sys/time.h>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include "mult.h"

typedef double number;

double get_time()
{
  struct timeval timeval_time;
  gettimeofday(&timeval_time,NULL);
  return (double)timeval_time.tv_sec + (double)timeval_time.tv_usec*1e-6;
}



int main(int argc, char *argv[])
{


  // define problem domain	
  const int dim = 3;
  const int nlocdofs = 1<<dim;
  printf("nlocdofs: %d \n",nlocdofs);

  int Ninterv = 1<<8;
      Ninterv = Ninterv - 6;
  printf("Ninterv: %d \n",Ninterv);
  const int N = 1+Ninterv;
  printf("N: %d \n",N);

  const int Nelems = Ninterv * (dim>1?Ninterv:1) * (dim>2?Ninterv:1);
  printf("Nelems: %d \n",Nelems);
  const int Ndofs = N * (dim>1?N:1) * (dim>2?N:1);
  printf("Ndofs: %d \n",Ndofs);

  // define solver's (pcg) convergence criteria
  int iters   = 0;
  int N_iters = 300;
  number tol  = 1e-5;
  
  // define some real matrices 
  number *u = new number[Ndofs];
  number *s = new number[Ndofs];
  number *f = new number[Ndofs];
  number *r = new number[Ndofs];
  number *p = new number[Ndofs];
  number *a = new number[Ndofs];
  number *M = new number[Ndofs];
  number *result = new number[Ndofs];
  number *vector = new number[Ndofs];

  // define local variables
  number *A   = new number[Nelems];
  number* Ke  = new number[nlocdofs*nlocdofs];
  number* prc = new number[nlocdofs*nlocdofs];

  // define misc iterative vars
  number *rho   = new number[N_iters];
  number *gamma = new number[N_iters];
  number *phi   = new number[N_iters];
  number *beta  = new number[N_iters];
  number *alpha = new number[N_iters];
  // initialize rho0,gamma0,etc
  rho[0]   = 0.0;
  gamma[0] = 0.0;
  phi[0]   = 0.0;
  beta[0]  = 0.0;
  alpha[0] = 0.0;

  // initialize system LHS & RHS (f[N])
  for(int i=0; i<Ndofs; i++){
	f[i] = 1e-5;
        s[i] = 0.0;
        u[i] = 1.0; //this value will change to u=0.0 after preconditioner
        a[i] = 0.0; //this value will change to u=0.0 after preconditioner
        s[i] = 0.0; //this value will change to u=0.0 after preconditioner
        p[i] = 0.0; //this value will change to u=0.0 after preconditioner
        r[i] = 0.0; //this value will change to u=0.0 after preconditioner
        M[i] = 0.0; //this value will change to u=0.0 after preconditioner
	if (i < 250) f[i] = 0.0;
  }
  
  // read local stiffness matrix Ke from "Ke.text" file
  std::fstream myfile("Ke.txt", std::ios_base::in);
  for (int i = 0; i <nlocdofs; i++){
      for (int j = 0; j <nlocdofs; j++){
	  int idx = (i*nlocdofs) + j;
        myfile >> Ke[idx];
      }
  }
	myfile.close();

  // simple Jacobi preconditioner [inverse of Kediag]
  for (int k = 0; k <nlocdofs*nlocdofs; k++){
      if ( (k % 9) == 0)
      prc[k] = Ke[k];	      
  }   

  // update jacobi preconditioner to global dofs
  SpMV(M,u,result,vector,prc,A,Ninterv);

   
  // Initialize vector u = u0 = 1/2^10=1/1024
  for (int i = 0; i <Ndofs; i++){
       u[i] = 1e-6;
       if ( (i % 9) == 0) M[i] = 1.0/M[i];
       //printf(" M %f\n",M[i]);
  }   

  // compute s[Ndofs]
  SpMV(s,u,result,vector,Ke,A,Ninterv);

  // update dot products
  for (int i = 0; i <Ndofs; i++){
     r[i] = f[i] - s[i];

     rho[0]   += M[i]*r[i]*r[i];
     gamma[0] += f[i]*f[i];

     p[i] = r[i];
  }   

  
/*** PCG Iterations Start Here ***/
double t = get_time();
while (sqrt(rho[iters]) > tol*sqrt(gamma[iters])){

  iters ++;

  // update p
  SpMV(a,p,result,vector,Ke,A,Ninterv);

  phi[iters] = 0.0;

  for (int i=0; i <Ndofs; i++) {
      phi[iters] += a[i]*p[i];
  }

  alpha[iters] = rho[iters - 1]/phi[iters];

  
  for (int i=0; i <Ndofs; i++) {
     u[i] += alpha[iters]*p[i];
     r[i] -= alpha[iters]*a[i];

     rho[iters] += M[i]*r[i]*r[i];
  }

  beta[iters] = rho[iters]/rho[iters-1];

  for (int i=0; i <Ndofs; i++) {
      p[i] = r[i] + beta[iters]*p[i];
  }

     printf("T  Residual %f    %f\n",u[10000000], sqrt(rho[iters]));

  if (iters == N_iters) {
     break;
  }
}
printf("Elaped Time: %g s\n",get_time()-t);




  delete[] u;
  delete[] s;
  delete[] f;
  delete[] r;
  delete[] p;
  delete[] a;
  delete[] M;
  delete[] A;
  delete[] result;
  delete[] vector;
  delete[] prc;
  delete[] rho;
  delete[] gamma;
  delete[] phi;
  delete[] beta;
  delete[] alpha;


  return 0;
}
