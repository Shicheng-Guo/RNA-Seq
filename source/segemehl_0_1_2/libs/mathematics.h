#ifndef _MATHEMATICS_H_
#define _MATHEMATICS_H_

#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"
#include <math.h>


#define RANDUNIT                        ((rand()/(RAND_MAX + 1.0)))
#define RANDINT(MAX)                    ((Uint)round(((double)MAX) * (rand()/((double)RAND_MAX + 1.0)))) 

#define MATRIX2D(X,NCOL,I,J) 			X[(I)*(NCOL)+(J)]
#define MATRIX3D(X,DIM_M,DIM_N,I,J,K) 	X[((I)*(DIM_M)+(J))*DIM_N+(K)]
#define MATRIX4D(X,DIM_M,DIM_N,DIM_Z,I,J,K,L) 	X[(((I)*(DIM_M)+(J))*DIM_N+(K))*DIM_Z+(L)]
#define VECTOR(X,I) 					((X)->elements[I])


#define INITMATRIX(X,SIZE,TYPESIZE) 	initArray(X,SIZE,TYPESIZE)
#define INITMATRIX2D(X,M,N,TYPESIZE) 	initArray(X, (M)*(N), TYPESIZE)
#define INITMATRIX3D(X,M,N,L,TYPESIZE) 	initArray(X, (M)*(N)*(L), TYPESIZE)
#define INITVECTOR(V) 					(V)->elements=NULL;\
										(V)->length=0

#define RESIZEVEC(V,N)					initArray(V,N,sizeof(vectorelem))
#define LENGTHVEC(V) 					((V)->length)
#define SWEEPVEC(V)						{int m;\
										for (m=0;m<LENGTHVEC(V);m++)\
										VECTOR(V,mi)=0;}
#define APPENDVEC(S,V,E) 				appendvector(S,V,E)
#define SWAPVEC(I,J,V) 	  				{ vectorelem msv; \
						  				msv = (V)->elements[I];\
						 				(V)->elements[I] = (V)->elements[J];\
						  				(V)->elements[J] = msv; }
#define DESTRUCTVEC(S,V)                destructVector((S),(V));
#define SWAPUINT(X,A,B)					{Uint msu ;\
  										 msu= X[A]; \
										 X[A] = X[B]; \
										 X[B]=msu; }

#define REVERSEVEC(A,B,V) 				{Uint mi;\
										for(mi=0; mi<(B-A); mi++) {\
										SWAPVEC(A+mi,B-mi,V);\
										}}
#define EMPTYVEC(V)						(V->elements==NULL)
#define MAX(A,B)						(((A) >= (B)) ? (A) : (B))
#define MAX3(A, B, C)					MAX(MAX((A),(B)),(C))
#define MAX4(A, B, C, D)				MAX(MAX(MAX((A),(B)),(C)),(D))
#define MAX5(A, B, C, D, E)				MAX(MAX(MAX(MAX((A),(B)),(C)),(D)),(E))

#define MIN(A,B)						(((A) <= (B)) ? (A) : (B))
#define MAX_DOUBLE						1e200
#define ABS(X)							((X)>=0?(X):-(X))
#define OVERLAP(A,B,C,D)                (((A) >= (C) && (B) <= (D)) || \
                                        ((A) <= (D) && (B) >= (D))  || \
                                        ((B) >= (C) && (A) <= (C)))            
#define CLOSEDINTERVAL(X,A,B)           ((A) >= (X) && (X) <= (B))
#define OPENINTERVAL(X,A,B)             ((A) > (X) && (X) < (B))
#define ISELEM(X,A,B)                   ((A) <= (X) && (X) <= (B))


#ifndef FLT_EPSILON
	#define FLT_EPSILON     1.192092896e-06 
	#define DBL_EPSILON     2.2204460492503131e-016
#endif

#ifndef M_PI
    #define M_PI 3.141592653589793238462643
#endif

#ifndef M_SQRT2
    #define M_SQRT2 1.4142135623730950488016887
#endif

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif

#ifndef VECTORELEM
typedef Uint vectorelem;
#endif

typedef struct {
	vectorelem *elements;
	Lint length;
} vector_t;

typedef struct{
  Uint noofbreaks;
  Uint *breaks;
  double RSS;
  double LL;
  double BIC;
} breakpoints_t;

void *initArray(void *, int, size_t);
void dumpMatrix_int(int *, int, int);
void dumpMatrix_Uint(Uint *, Uint, Uint);
void dumpMatrix_flt(float *, int, int);
void dumpMatrix_dbl(double *, Uint, Uint);
void dumpMatrix3D_int(int *, int, int, int);
double *transpose(void* space, double *, Uint, Uint);
void appendvector(void *, vector_t *, vectorelem);
void dumpVector(vector_t *);
void destructVector(void *, vector_t *);
Uint uarraymax(Uint *, Uint);
int arraymax(int *, int);
int gcd(int, int);
double power(double, int);
double uniroot(double, double, double (*f)(double, void*), double, void*);
double BLAST_Expm1(double x);
Uint fak(Uint);
double* coldel (void *, double *, Uint, Uint, Uint);
double* rowdel (void *, double *, Uint, Uint, Uint);
Uint minvecdist(void *space, vector_t *vec, Uint i, Uint j);
int* intrev(int *n, Uint len);
double shannonentropy(void *space, char *seq, Uint len, Uint asize, Uint *encodetab);
double log10(double x);
double logadd(double a, double b);
double log10add(double a, double b);
double var_int (int *x, Uint n);
double poisson(double lambda, double x);
double logpoisson(double lambda, double x);
double multivarnorm (double *pt, double *mu, double *sd, Uint n);
double bivarnorm(double x, double y, double mu1, double mu2, double* cv);
double univarnormcdf (double x, double mu, double sd);
double randunivarnorm(double mu, double sd);
double gmm(void *space, double *pt, Uint m, Uint n, double *mu, double *sd, double *w, Uint g, Uint maxiter);
double var (double *x, Uint n);
void normalize (double *a, Uint n);
breakpoints_t * bl_RSSmatrix (void *space, double *x, Uint n, Uint h, Uint noofbreaks);
double* bl_RSS (void *space, double *x, Uint n, Uint u, Uint v);


#endif
