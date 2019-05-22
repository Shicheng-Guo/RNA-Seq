
/*
 * mathematics.c
 * implemtation of various mathematical functions
 *
 * @author Steve Hoffmann
 * @date Wed 22 Nov 2006
 *
 *  SVN
 *  Revision of last commit: $Rev: 54 $
 *  Author: $Author: steve $
 *  Date: $Date: 2008-09-10 22:13:30 +0200 (Wed, 10 Sep 2008) $
 *
 *  Id: $Id: mathematics.c 54 2008-09-10 20:13:30Z steve $
 *  Url: $URL: http://www.bioinf.uni-leipzig.de/svn/segemehl/segemehl/trunk/libs/mathematics.c $
 *  
 */
#include "mathematics.h"
#include <string.h>
#include <limits.h>
#include <math.h>


int* intrev(int* n, Uint len){
  int end = len-1;
  int start = 0;

  while (start<end) {
    n[start] ^= n[end];
    n[end] ^= n[start];
    n[start] ^= n[end];
    start++;
    end--;
  }
  return n;
}

 void *initArray(void *space, int size, size_t datatype) {
	void *ptr=NULL;

	/*dirty trick: sizeof(char) == 1*/
	ptr = ALLOCMEMORY(space, ptr, char, size*datatype);
	return ptr;
 }


void appendvector(void *space, vector_t *v, vectorelem elem) { 

  	 v->elements = (vectorelem*) ALLOCMEMORY(space, v->elements, vectorelem, (v->length+1));
	 v->elements[v->length]=elem;
	 v->length++;
}




/*--------------------------------- mindist ----------------------------------
 *    
 * @brief expects a sorted vector to find the minimum distance between 
 * vec[i] and vec[j]
 * @author Steve Hoffmann 
 *   
 */

Uint
minvecdist(void *space, vector_t *vec, Uint i, Uint j) {
  Uint k, 
       size_i,
       size_j;
  int range,
      dist = INT_MAX,
      l;
  vectorelem *e_i,
             *e_j;

  size_j = LENGTHVEC(&vec[j]);
  size_i = LENGTHVEC(&vec[i]);

  if (size_i == 0 || size_j == 0) 
    return 0;

  e_j = &vec[j].elements[0];
  for(k=0; k < size_j ; k++, e_j++) {
    e_i = &vec[i].elements[0];        
    for(l=0; l < size_i; l++, e_i++) {
      range = abs((int)*e_j - (int)*e_i);
      if (range < dist) {
        dist = range;
      }
    }
  }

  return dist;
}


Uint
minvecdist2(void *space, vector_t *vec1, vector_t *vec2, Uint *which) {
  Uint k, 
       size_i,
       size_j;
  int range,
      dist = INT_MAX,
      l;
  vectorelem *e_i,
             *e_j;

  size_j = LENGTHVEC(vec2);
  size_i = LENGTHVEC(vec1);

  if (size_i == 0 || size_j == 0) 
    return 0;

  e_j = &vec2->elements[0];
  for(k=0; k < size_j ; k++, e_j++) {
    e_i = &vec1->elements[0];        
    for(l=0; l < size_i; l++, e_i++) {
      range = abs((int)*e_j - (int)*e_i);
      if (range < dist) {
        dist = range;
        *which = l;
      }
    }
  }

  return dist;
}



void dumpMatrix_int(int *M, int m, int n) {
	int i,j;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			printf("%d ", MATRIX2D(M,n,i,j));
		}
			printf("\n");
	}
 }

Uint uarraymax(Uint *arr, Uint l) {
	Uint i;
	Uint max =0;
	
  	for(i=0; i < l; i++) {
		if (arr[i]>arr[max]) max=i;
	}
	
	return max;
}



int arraymax(int *arr, int l) {
	int i;
	int max =0;
	
  	for(i=0; i < l; i++) {
		if (arr[i]>arr[max]) max=i;
	}
	
	return max;
}

void dumpMatrix_Uint(Uint *M, Uint m, Uint n) {
	Uint i,j;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			printf("%d ", MATRIX2D(M,n,i,j));
		}
			printf("\n");
	}
 }


void dumpMatrix_dbl(double *M, Uint m, Uint n) {
	Uint i,j;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			printf("%f ", MATRIX2D(M,n,i,j));
		}
			printf("\n");
	}
 }


 void dumpMatrix3D_int(int *M, int m, int n, int l) {
	int i,j,k;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			for (k=0; k < l; k++) {
				printf("%d ", MATRIX3D(M,n,l,i,j,k));
			}
		printf(";");
		}
	printf("\n");
	}
 }

void dumpVector(vector_t *v) {

	int i;
	for (i=0; i < v->length; i++) {
		printf("%d ", v->elements[i]);
	}

	printf("\n");
}


void destructVector(void *space, vector_t *v) {
	
    if (v!=NULL) {
    	if (v->elements) FREEMEMORY(space, v->elements);
		FREEMEMORY(space, v);
	}
}

void reverseVector(Uint a, Uint b, vector_t *v) {
	Uint i;
	
	for (i=0; i < (b-a); i++) {
		SWAPVEC(a+i,b-i,v);
	}
}

int nextPermutation(vector_t *v) {
	Uint i,j; 
	vectorelem *e=v->elements;

	for (i=(v->length)-1; i > 0; i--)
		if(e[i-1]<=e[i]) break;
	
	if (i==0) return 0;

	for (j=i+1; j < (Uint) v->length; j++ )
		if(e[i-1]>=e[j]) break;
	
	SWAPVEC(i-1, j-1, v);
	REVERSEVEC(i, (v->length)-1, v);

	return 1;
}

/*----------------------------------- norm -----------------------------------
 *    
 * @brief normalize
 * @author Steve Hoffmann 
 *   
 */
 
void
normalize (double *a, Uint n)
{
  Uint i;
  double sum = 0;

  for(i=0; i < n; i++) {
    sum += a[i];
  }
	  
  for(i=0; i < n; i++) {
    a[i] /= sum;
  }
 
  return ;
}



/*----------------------------------- gcd ------------------------------------
 *    
 * calculate the greatest common divisor of two integer values
 * 
 */
 
int
gcd (int a, int b)
{
    int val;

	b = abs(b);
	
	if (b > a)
	  val=a, a=b, b=val;

	while (b != 0) {
		val = a%b;
		a = b;
		b = val;
	}
	
	return a;
}


/*---------------------------------- power -----------------------------------
 *    
 * the power may be with you! 
 * 
 */
 
double
power (double x, int n)
{
  	double y = 1.;

	if(n==0)
	  return 1;
	if(x==0) {
		if(n < 0) {
			return MAX_DOUBLE;
		}
		return 0;
	}

	if (n < 0) {
		x = 1./x;
		n = -n;
	}

	while(n > 0) {
		if (n & 1) {
			y *= x;
		}
		n /= 2;
		x *= x;
	}

	return y;
}




/*----------------------------------- fak ------------------------------------
 *    
 * @brief get the factorial (works only for n <=10!!)
 * @author Steve Hoffmann 
 *   
 */
 

Uint fak(Uint n) {
  Uint i,x=n;
  
  for(i=x-1; i > 0; i--) {
  	x *= i;
  }

  return x;
}


/*--------------------------------- uniroot ----------------------------------
 *    
 * getting the zero-root of a given function
 * 
 * according to G. Forsythe, M. Malcom et al.
 * Computer methods for mathematical computations, 1980
 * 
 */
 
double
uniroot (double start, double end, double (*f)(double, void*), 
    double tolx, void* info)
{	
  	double a, b, c;
	double fa, fb, fc;
	double prev;
	double currenttol;
	double p, q, new_step;
	double cb, t1, t2;

	a = start; b= end; fa = (*f)(a,info); fb=(*f)(b,info);
	c = a; fc = fa;
	
	if ((fa > (double) 0 && fb > (double) 0) 
        || (fa < (double)0 && fb < (double)0)) {
		printf("mooep!\n");	
	  /*return 0;*/
	} 
	
	while(1) {

	  	prev = b-a;
		
		if (fabs(fc) < fabs(fb)) {
			a=b; b=c; c=a;
			fa=fb; fb=fc; fc=fa;
		}
		currenttol = 2 * FLT_EPSILON * fabs(b) + tolx/2;
		new_step = (c-b)/2;
		if (fabs(new_step) <= currenttol || fb == (double)0) {
			return b;
		}
		
		if ( fabs(prev) >= currenttol && fabs(fa) > fabs(fb) ) {
			cb = c-b;
			if(a==c) {
				t1 = fb/fa;
				p = cb*t1;
				q = 1.0 - t1;
			} else {
				q = fa/fc;
				t1 = fb/fc;
				t2 = fb/fa;
				p = t2 * ( cb * q * (q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);	
			}
			if ( p > (double)0) {
				q = -q;
			} else {
				p = -p;
			}

			if(p < (0.75 * cb * q - fabs(currenttol*q)/2) 
				&& p < fabs(prev * q/2) ) {
				new_step = p/q;
			}
		}

		if (fabs(new_step) < currenttol ) {
			if(new_step > (double)0) {
				new_step = currenttol;
			} else {
				new_step = -currenttol;
			}
		}
		
		a=b; fa=fb;
		b+= new_step;
		fb = (*f)(b,info);
		if( (fb>0 && fc>0) || (fb < 0 && fc < 0) ) {
			c=a; fc=fa;
		}	
	}
	
	return 0;
}


/*---------------------------------- coldel ----------------------------------
 *    
 * @brief delete column for matrix
 * @author Steve Hoffmann 
 *   
 */
 
double*
coldel (void *space, double *a, Uint m, Uint n, Uint d) {
	
	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=0;

  t = (double*) INITMATRIX2D(space, m, (n-1), sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) { 
	  j++; k=0; l=0;
	} 	
	if(k++ != d) {
	  MATRIX2D(t, n-1, j, l++) = a[i];
	}
  }
	
  FREEMEMORY(space, a);
  return t;
}


/*---------------------------------- rowdel ----------------------------------
 *    
 * @brief delete row from matrix
 * @author Steve Hoffmann 
 *   
 */
 
double*
rowdel (void *space, double *a, Uint m, Uint n, Uint d) {
	
	double *t;
	Uint	i,
			j=-1,
			k=0,
			l=-1;

  t = (double*) INITMATRIX2D(space, (n-1), m, sizeof(double));

  for(i=0; i < m*n; i++) {
	if(i % n == 0) { 
	  j++; k=0;
	  l = (j != d) ? l+1 : l;
	} 	
	if(j != d) {
	  MATRIX2D(t, n, l, k++) = a[i];
	}
  }

  FREEMEMORY(space, a);
  return t;
}


/*---------------------------------- xprod -----------------------------------
 *    
 * @brief calculate the cross product of two vectors
 * @author Steve Hoffmann 
 *   
 */
 
double*
xprod (void *space, double* x, Uint m, double *y, Uint n) {
	double *p;
	Uint 	i,	
			j;

	p = (double*) INITMATRIX2D(space, m, n, sizeof(double));

	for (i=0; i < m; i++) {
		for(j=0; j < n; j++) {
			MATRIX2D(p, n, i, j) = x[i]*y[i];
		}
	}
	return p;
}


/*-------------------------------- transpose ---------------------------------
 *    
 * @brief transpose a matrix $a$ of dimensions $m x n$
 * @author Steve Hoffmann 
 *   
 */

double*
transpose (void* space, double *a, Uint m, Uint n) {
  double *t;
  Uint	i,
        j=-1,
        k=0;

  t = (double*) INITMATRIX2D(space, n, m, sizeof(double));

  for(i=0; i < m*n; i++) {
    if(i % n == 0) { j++; k=0;} 	
    MATRIX2D(t, m, k, j) = a[i];
    k++;
  }

  FREEMEMORY(space, a);
  return t;
}


/*--------------------------------- simpson ----------------------------------
 *    
 * @brief implementation of the simpson algorithm to determine the integral
 * of $f(x)$ in the interval [$a$,$b$]. sdiv denotes number of subdivisions
 * @author Steve Hoffmann 
 *   
 */

  double
simpson( double a, double b, int sdiv, 
    double (*f) (double, void*),
    void* info) 
{

  double 	k,
            sum1=0,
            sum2=0;
  int 	i;

  k = ((double) b-a)/((double)2*sdiv);

  for (i=1; i < sdiv; i++) {
    sum1+=f(a + k*2*i, info);
    sum2+=f(a + k*(2*i-1), info);
  }

  sum2+=f(a + k*(2*i-1), info);

  return ((double)(k/3)*(f(a, info)+f(b, info)+2*sum1+4*sum2));
}


/*-------------------------------- simpson1D ---------------------------------
 *    
 * @brief helper function for simpson2D
 * @author Steve Hoffmann 
 *   
 */

  double
simpson1D(double x, int sdiv, 
    double (*f) (double, double, void*),
    double (*c) (double, void*),
    double (*d) (double, void*),
    void* info) 
{

  double 	k,
            sum1=0,
            sum2=0,
            ca,
            da;
  int 	    i;

  ca = c(x, info);
  da = d(x, info);

  k = ((double) da-ca)/((double)2*sdiv);

  for (i=1; i < sdiv; i++) {
    sum1+=f(x, ca + k*2*i, info);
    sum2+=f(x, ca + k*(2*i-1), info);
  }

  sum2+=f(x, ca + k*(2*i-1), info);


  return ((double)(k/3)*(f(x, ca, info)+f(x, da, info)+2*sum1+4*sum2));
}


/*-------------------------------- simpson2D ---------------------------------
 *    
 * @brief calculates the 2-dim integral of function $f$ given the interval
 * [$a$,$b$] in the first and [$c(x)$,$d(x)$] in the second dimension
 * sdiv, sdiv2 denote the subdivisions in the first and second dimension
 * @author Steve Hoffmann 
 *   
 */

double 
simpson2D(double a, double b, int sdiv, int sdiv2, 
    double (*f) (double, double, void*), 
    double (*c) (double, void*),
    double (*d) (double, void*),
    void *info) {

  double 	h,	
  sum1=0,
  sum2=0;
  int		i;

  h = ((double)b-a)/((double)2*sdiv);

  for (i=1; i < sdiv; i++) {
    sum1 += simpson1D((a+h*2*(i)), sdiv2, f, c, d, info);
    sum2 += simpson1D((a+h*(2*i-1)), sdiv2, f, c, d, info);
  }

  sum2 += simpson1D((a+h*(2*i-1)), sdiv2, f, c, d, info);

  return ((double)(h/3) * (simpson1D(a, sdiv2, f, c, d, info) + 
        simpson1D(b, sdiv2, f, c, d, info) + 2*sum1 + 4*sum2 ));
}



/*--------------------------------- myMinor ----------------------------------
 *    
 * @brief helper function for the laplacian algorithm used in det()
 * @author Steve Hoffmann 
 *   
 */
double* 
myMinor(void *space, double* M, Uint m, Uint n, Uint i, Uint j) {

  double *t;

  t = (double*) ALLOCMEMORY(space, NULL, double, m*n);
  memmove(t, M, sizeof(double)*(m*n));

  t = rowdel(NULL, t, m, n, i);
  t = coldel(NULL, t, m-1, n, j);  

  return t;
}



/*----------------------------------- det ------------------------------------
 *    
 * @brief calculates the determinant of a square matrix m of size n x n using 
 * Laplacian algorithm (recursive implementation)
 * @author Steve Hoffmann 
 *   
 */

double
det(void *space, double *m, int n) {
  double sum=0,
         *t=NULL;
  int		j;

  if (n==1) {	
    return MATRIX2D(m, n, 0, 0);
  }

  for(j=0; j < n; j++) {

    t = myMinor(space, m, n, n, 0, j);
    sum += pow(-1.0, (j+2))*MATRIX2D(m,n,0,j)*det(space, t, n-1);
  }

  FREEMEMORY(space, t);
  return sum;
}


/*----------------------------------- add ------------------------------------
 *    
 * @brief componentwise addition of a to a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double*
add(double *x, Uint m, double a) {
  Uint i;

  for(i=0; i < m; i++) {
    x[i] += a;
  }
  return x;
}


/*----------------------------------- mean -----------------------------------
 *    
 * @brief calculate the arithmetic mean for a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double
mean (double *x, Uint m) {
  Uint i;
  double sum=0;

  for (i=0; i < m; i++) {
    sum += x[i];
  }

  return sum /= m; 
}

/*----------------------------------- mean -----------------------------------
 *    
 * @brief calculate the arithmetic mean for a vector of length m
 * @author Steve Hoffmann 
 *   
 */

double
mean_int (int *x, Uint m) {
  Uint i;
  double sum=0;

  for (i=0; i < m; i++) {
    sum += x[i];
  }

  return sum /= m; 
}



/*---------------------------------- scalar ----------------------------------
 *    
 * @brief calculate the scalar product of two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double
scalar (double* x, double *y, Uint m) {
  double  p=0;
  Uint 	i;

  for (i=0; i < m; i++) {
    p += x[i]*y[i];
  }
  return p;
}


/*----------------------------------- cov ------------------------------------
 *    
 * @brief get the covariance matrix (2x2) for two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double*
cov (void *space, double *x, double *y, Uint m) {
  double *c,
         xm,
         ym;

  c = (double*) INITMATRIX2D(space, 2, 2, sizeof(double));
  xm = mean(x, m);
  ym = mean(y, m);

  /*center*/
  add(x, m, (-1)*xm);
  add(y, m, (-1)*ym);

  MATRIX2D(c, 2, 0, 0) = (double) scalar(x,x,m)/(m-1);
  MATRIX2D(c, 2, 0, 1) = MATRIX2D(c, 2, 1, 0) = (double) scalar(x,y,m)/(m-1);
  MATRIX2D(c, 2, 1, 1) = (double) scalar(y,y,m)/(m-1);

  return c;
}




/*----------------------------------- var ------------------------------------
 *    
 * @brief get the sample variance
 * @author Steve Hoffmann 
 *   
 */
 
double
samplevar (double *x, double *p, double n)
{   
    int i;
    double m, r, sum=0;

    m=mean(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r)*p[i];
    }

	return sum/n;
}


/*----------------------------------- var ------------------------------------
 *    
 * @brief get the variance
 * @author Steve Hoffmann 
 *   
 */
 
double
var_int (int *x, Uint n)
{   
    int i;
    double m, r, sum=0;

    m=mean_int(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r);
    }

	return sum/n;
}

/*--------------------------------- poisson ----------------------------------
 *    
 * @brief the <.><< distribution
 * @author Steve Hoffmann 
 *   
 */
 
double
poisson(double lambda, double x) {
  assert(x >= 0);
  return (pow(lambda,x)/tgamma(x+1))*exp(-lambda);
}


/*-------------------------------- logpoisson --------------------------------
 *    
 * @brief poisson in log space
 * @author Steve Hoffmann 
 *   
 */
 
double
logpoisson (double loglambda, double logx)
{
  //assert(logx >= 0);
  return (exp(logx)*loglambda) - log(tgamma(exp(logx)+1)) - exp(loglambda);   
	
}

/*----------------------------------- var ------------------------------------
 *    
 * @brief get the variance
 * @author Steve Hoffmann 
 *   
 */
 
double
var (double *x, Uint n)
{   
    int i;
    double m, r, sum=0;

    m=mean(x, n);
    for (i=0; i < n; i++) {
      r = x[i]-m;
      sum += (r*r);
    }

	return sum/n;
}


/*---------------------------------- stddev ----------------------------------
 *    
 * @brief get the standard deviation
 * @author Steve Hoffmann 
 *   
 */
 
double
stddev (double *x, double n)
{
	
    return sqrt(var(x, n));
}


/*----------------------------------- rho ------------------------------------
 *    
 * @brief calculate correlation $\rho$ for two vectors of length m
 * @author Steve Hoffmann 
 *   
 */

double
rho (void *space, double *x, double *y, Uint m) {
  double *cv;

  cv = cov(space, x, y, m); 
  return (MATRIX2D(cv, 2, 0, 1)/sqrt(MATRIX2D(cv, 2, 0, 0)*MATRIX2D(cv, 2, 1, 1)));
}

/*-------------------------------- univarnorm --------------------------------
 *    
 * @brief pdf gaussian
 * @author Steve Hoffmann 
 *   
 */
 
double
univarnorm (double x, double mu, double sd)
{
    double d = (x-mu);
    double sdsq = sd * sd;

	return exp(-0.5*(d*d/sdsq))/sqrt(2*M_PI*sdsq);
}


/*------------------------------ univarnormcdf -------------------------------
 *    
 * @brief cdf gaussian
 * @author Steve Hoffmann 
 *   
 */
 
double
univarnormcdf (double x, double mu, double sd)
{
	return 0.5 * (1+erf((x-mu)/(sd*M_SQRT2)));
}


/*------------------------------ randunivarnorm ------------------------------
 *    
 * @brief algorithm adapted from Dr. Everett (Skip) Carter, Jr.
 * @author Steve Hoffmann 
 *   
 */
 
double
randunivarnorm (double mu, double sd)
{
  double x1, x2, w, y1, y2;
  
  do{
    x1 = 2.0 * (((double)rand())/((double)RAND_MAX)) - 1.0;
    x2 = 2.0 * (((double)rand())/((double)RAND_MAX)) - 1.0;
    w = x1 * x1 + x2 * x2; 
  } while (w >= 1.0 || w == .0);

  w = sqrt((-2.0*log(w))/w);
  y1 = x1 *w;
  y2 = x2 *w;

  return y2*sd+mu;
}

/*-------------------------------- bivarcond ---------------------------------
 *    
 * @brief conditional bivar. norm. distrib. f(y|x) given location parameter
 * $mu1$, $mu2$ and covariance matrix $cv$ of size (2x2)
 * @author Steve Hoffmann 
 *   
 */

double
bivarcond(double x, double y, double mu1, double mu2, double *cv) {
  double rho,
         s1,
         s1sq,
         s2,
         s2sq,
         m,
         e;

  s1sq = MATRIX2D(cv, 2, 0, 0);
  s2sq = MATRIX2D(cv, 2, 1, 1);
  s1 = sqrt(s1sq);
  s2 = sqrt(s2sq);
  rho = MATRIX2D(cv, 2, 0, 1)/sqrt(s1sq*s2sq);

  m  = 1/sqrt((2*M_PI*s2sq*(1-(rho*rho))));
  e = (y-mu2-rho*(s2/s1)*(x-mu1));
  e *= e;
  e /= s2sq*(1-(rho*rho));

  return(m*exp(-0.5*e));
}


/*-------------------------------- bivarnorm ---------------------------------
 *    
 * @brief bivariate normal distribution f(x,y) given location parameter
 * $mu1$, $mu2$ and covariance matrix $cv$ of size (2x2)
 * @author Steve Hoffmann 
 *   
 */


double
bivarnorm(double x, double y, double mu1, double mu2, double* cv) {
  double rho,
         s1,
         s1sq,
         s2,
         s2sq,
         m,
         e1,
         e2;

  s1sq = MATRIX2D(cv, 2, 0, 0);
  s2sq = MATRIX2D(cv, 2, 1, 1);
  s1 = sqrt(s1sq);
  s2 = sqrt(s2sq);
  rho = MATRIX2D(cv, 2, 0, 1)/sqrt(s1sq*s2sq);

  m = 1/(2*M_PI*s1*s2*sqrt(1-(rho*rho)));

  e1 = (-1)/(2*(1-(rho*rho)));
  e2 = ((x-mu1)*(x-mu1))/s1sq  
    - (2*rho*(x-mu1)*(y-mu2))/(s1*s2)
        + ((y-mu2)*(y-mu2))/s2sq; 

  return m*exp(e1*e2);
}

/*------------------------------ multivarnorm -------------------------------
 *    
 * @brief n-dimensional gaussian probability density function wo correlation term
 * i.e. orthogonal variation!
 * @author Steve Hoffmann 
 *   
 */
 
double
multivarnorm (double *pt, double *mu, double *sd, Uint n)
{
  Uint i;
  double det = 1; 
  double exponent = 0;
  double perturb = 0;

  for(i=0; i < n; i++) {
    exponent += pow((pt[i]-mu[i]),2) * (1/pow(MAX(sd[i],perturb),2));
    det *= pow(sd[i],2);
  }

  //fprintf(stderr, "exponent: %f, det: %f, exp(expo):%f, norm:%f\n", exponent, det,
  //    exp(exponent*-0.5), sqrt(pow((2*M_PI),n)*det));

  return exp(exponent*-0.5)/sqrt(pow((2*M_PI),n)*det);	
}


/*this is ncbi intellectual property*/

double BLAST_Expm1(double x)
{
  double	absx = ABS(x);

  if (absx > .33)
	return exp(x) - 1.;

  if (absx < 1.e-16)
	return x;

  return x * (1. + x *
	  (1./2. + x * 
	   (1./6. + x *
		(1./24. + x * 
		 (1./120. + x *
		  (1./720. + x * 
		   (1./5040. + x *
			(1./40320. + x * 
			 (1./362880. + x *
			  (1./3628800. + x * 
			   (1./39916800. + x *
				(1./479001600. + 
				 x/6227020800.))))))))))));
}


double
log10(double x) {
  return (log(x)/log(10));
}

double
log2(double x) {
  return (log(x)/log(2));
}

double 
log10add(double a, double b) {
  double max, min;

  if(a > b) {
    if (b == log10(0)) return a;
    else {
      max=a;
      min=b;
    }
  } else {
    if (a == log10(0)) return b;
    else {
      max=b;
      min=a;
    }
  }
  return max + log10(1+pow(10.,min-max));  
}

double 
logadd(double a, double b) {
  double max, min;

  if(a > b) {
    if (b == log(0)) return a;
    else {
      max=a;
      min=b;
    }
  } else {
    if (a == log(0)) return b;
    else {
      max=b;
      min=a;
    }
  }
  return max + log(1+exp(min-max));  
}

/*-------------------------------- seqentropy --------------------------------
 *    
 * @brief zero order sequence entropy
 * @author Steve Hoffmann 
 *   
 */
 

double
shannonentropy(void *space, char *seq, Uint len, Uint asize, Uint *encodetab) {
  //CHANGED: count should start with 0
  //--> possibly just use len instead of norm 
  //Uint i, norm=1;
  Uint i, norm=0;
  double *p, H=0;


  p = ALLOCMEMORY(space, NULL, double, asize);
  memset(p, 0, sizeof(double)*asize);

  for(i=0; i < len; i++) {
    p[encodetab[(Uint)seq[i]]]++;
    norm++;
  }

  for(i=0; i < asize; i++) {
    p[i]/=norm;
  }

  for(i=0; i < asize; i++) {
    if(p[i] > 0)
    H += p[i] * log2(p[i]); 
  }

  FREEMEMORY(space, p);
  return -1*H;
}

/*-------------------------------- smoothavg ---------------------------------
 *    
 * @brief smoothing average
 * @author Steve Hoffmann 
 *   
 */
 
double*
smoothavg (double *y, Uint n, Uint r)
{
  double *ys;
  Uint i, j;

  ys = ALLOCMEMORY(space, NULL, double, n);

  for(i=r; i < n; i++) {
    for(j=0; j < 2*r; j++) {
      ys[r+j] += y[i-r+j] / (2*r+1);
    }
  }

  return ys;
}

/*----------------------------------- gmm ------------------------------------
 *    
 * @brief gaussian mixture model for m n-dimensional data points and g gaussians
 * @author Steve Hoffmann 
 *   
 */
 
double
gmm(void *space, double *pt, Uint m, Uint n, 
    double *mu, double *sd, double *w, Uint g, Uint maxiter) {

  double *ms, *mu_, *sd_, no, dd, ll, oll, epsilon=0.00000000001;
  Uint i, j, k, l=0;

  ms = ALLOCMEMORY(space, NULL, double, m*g);
  mu_ = ALLOCMEMORY(space, NULL, double, n);
  sd_ = ALLOCMEMORY(space, NULL, double, n);
  memset(sd_, 0, sizeof(double)*n);
  memset(mu_, 0, sizeof(double)*n);
  memset(ms, 0, sizeof(double)*m*g);

  ll = 1;

  do {

    oll = ll;

    /*expectation step*/
    for(ll=0, i=0; i < m; i++) {
      for(j=0; j < g; j++) {
        MATRIX2D(ms, g, i, j) = 
          multivarnorm(&MATRIX2D(pt, n, i, 0), &MATRIX2D(mu, n, j, 0),
              &MATRIX2D(sd, n, j, 0), n);
        ll += log(w[j] * MATRIX2D(ms, g, i, j));
      }
      normalize(&MATRIX2D(ms, g, i, 0), g);
    }

    /*maximization step - weighted normalized mean*/
    for(j=0; j < g; j++) {
      memset(mu_, 0, sizeof(double)*n);

      for(i=0, no=0; i < m; i++) {
        no += MATRIX2D(ms, g, i, j);
        for(k=0; k < n; k++) {
          mu_[k] += MATRIX2D(pt, n, i, k) * MATRIX2D(ms, g, i, j);
        }
      }

      /*update mean*/   
      for(k=0; k < n; k++) { 
        mu_[k] /= no;
        MATRIX2D(mu, n, j, k) = mu_[k];
      }

      for(i=0; i < m; i++) {
        for(k=0; k < n; k++) {
          dd = MATRIX2D(pt, n, i, k) - mu_[k];
          sd_[k] += (dd * dd) * MATRIX2D(ms, g, i, j);
        }
      }

      /*update sd*/
      for(k=0; k < n; k++) { 
        sd_[k] /= no;
        MATRIX2D(sd, n, j, k) = sqrt(sd_[k]);
      }
    }
   

  } while(l++ < maxiter && fabs((double)ll-oll) > -1*ll*epsilon);

  FREEMEMORY(space, mu_);
  FREEMEMORY(space, sd_);
  FREEMEMORY(space, ms);

  return ll;
}



/*--------------------------------- bl_RSSi ----------------------------------
 *    
 * @brief calculation of a RSS for a vector of length n from u to v for intercept
 * @author Steve Hoffmann 
 *   
 */
 
double*
bl_RSS (void *space, double *x, Uint n, Uint u, Uint v)
{
  Uint i;
  double *cum, *y;
  
  assert(v>u);
  
  cum = ALLOCMEMORY(space, NULL, double, v-u+1);
  y = ALLOCMEMORY(space, NULL, double, v-u+1);

  cum[0] = x[u];
  y[0] = x[u] * sqrt(2);

  for(i=u+1; i < v; i++) {
    cum[i-u] = x[i] + cum[i-u-1]; 
    y[i-u] = x[i] - cum[i-u]/((double)i-u+1);
    y[i-u] *= sqrt(1.0+(1.0/((double)i-u)));
  }
      
  cum[0] = 0;

  for(i=1; i < v-u+1; i++) {
    cum[i] = cum[i-1] + y[i]*y[i];
  }

  //diagonal
  cum[0] = 0;

  FREEMEMORY(space, y);
  return cum;
}


/*------------------------------- bl_RSSmatrix -------------------------------
 *    
 * @brief compute the triangular RSS matrix 
 * @author Steve Hoffmann 
 *   
 */
 
breakpoints_t*
bl_RSSmatrix (void *space, double *x, Uint n, Uint h, Uint noofbreaks)
{
  Uint i, j, minarg, *POS, m, resc;  //m breaks
  breakpoints_t *breaks;
  double *M, *RSS, *temp, minval;
  double *y;

  /*calculation of rss matrix*/
  M = ALLOCMEMORY(space, NULL, double, (n-h+1)*n);
  memset(M, 0, sizeof(double)*((n-h+1)*n));

  for(i=0; i < n-h+1; i++) {
    y = bl_RSS(space, x, n, i, n);
    memmove(&MATRIX2D(M, n, i, i), y, sizeof(double)*(n-i));
    FREEMEMORY(space, y);
  }
  
  /*DP in the rss matrix*/
  RSS = ALLOCMEMORY(space, NULL, double, (noofbreaks+1)*(n+1));
  POS = ALLOCMEMORY(space, NULL, Uint, (noofbreaks+1)*(n+1));
  memset(RSS, 0, (noofbreaks+1)*(n+1)*sizeof(double));
  memset(POS, 0, (noofbreaks+1)*(n+1)*sizeof(Uint));

  for(i=0; i < n; i++) {
    MATRIX2D(RSS, noofbreaks+1, i, 1) = MATRIX2D(M, n, 0, i);
    MATRIX2D(POS, noofbreaks+1, i, 1) = i;
  }

  for(m=2; m <= noofbreaks; m++) { 
    for(i=m*h-1; i < n-h; i++) {
      temp = ALLOCMEMORY(space, NULL, double, (i-h)-((m-1)*h-1)+1);
      /*Recursion*/
      for(j=(m-1)*h-1; j < i-h+1; j++) {
        temp[j-(((m-1)*h)-1)] = MATRIX2D(RSS, noofbreaks+1, j, m-1) + MATRIX2D(M, n, j+1, i);
      }
      /*find the minimum*/
      minarg = 0;
      minval = temp[0];
      for(j=1; j < (i-h)-((m-1)*h-1)+1; j++) {
        if(temp[j] < minval) {
          minarg = j;
          minval = temp[j];
        }
      }
      /*register*/
      MATRIX2D(RSS, noofbreaks+1, i, m) = temp[minarg];
      MATRIX2D(POS, noofbreaks+1, i, m) = (m-1)*h-1+minarg;

      FREEMEMORY(space, temp);
    }
  }

  
  breaks = ALLOCMEMORY(space, NULL, breakpoints_t, noofbreaks+1);
  breaks[0].noofbreaks = 0;
  breaks[0].RSS = MATRIX2D(RSS, noofbreaks+1, n-1, 1);
  breaks[0].LL = -0.5 * n * (log(breaks[0].RSS) + 1 - log(n) + log(2*M_PI));
  breaks[0].BIC = -2*breaks[0].LL + 2*log(n);
  breaks[0].breaks = NULL;


  for(m=noofbreaks; m >= 1; m--) { 
    /*calling breaks and backtrace*/
    temp = ALLOCMEMORY(space, NULL, double, n);

    for(i=h-1; i < n-h; i++) {
      if (MATRIX2D(RSS, noofbreaks+1, i, m) == 0.0 || 
          MATRIX2D(M, n, i+1, n-1) == 0.0) {
        temp[i]= -1.0;
      } else {  
        temp[i] = MATRIX2D(RSS, noofbreaks+1, i, m) + MATRIX2D(M, n, i+1, n-1);
      }
    }

    minarg = 0;
    minval = temp[h-1];
    for(i=h-1; i < n-h; i++) {
      if((temp[i] > 0.0 && temp[i] < minval) || minval == -1.0) {
        minarg = i;
        minval = temp[i];
      }
    }

    breaks[m].noofbreaks = m; 
    breaks[m].breaks = ALLOCMEMORY(space, NULL, Uint, m);
    breaks[m].RSS = temp[minarg];
    breaks[m].LL =  -0.5 * n * (log(temp[minarg]) + 1 - log(n) + log(2*M_PI));
    breaks[m].BIC =  -2*breaks[m].LL + (2*(m+1))*log(n);
    breaks[m].breaks[0] = minarg;


    for(i=m, j=1; i >= 2; i--, j++) {
       breaks[m].breaks[j] = MATRIX2D(POS, noofbreaks+1, breaks[m].breaks[j-1], i);
    }
 
    //reverse
    for(i=0; i < m/2; i++) { 
      resc = breaks[m].breaks[m-i-1];
      breaks[m].breaks[m-i-1] = breaks[m].breaks[i];
      breaks[m].breaks[i] = resc;
    }
    FREEMEMORY(space, temp);
  }
   
  FREEMEMORY(space, RSS);
  FREEMEMORY(space, M);
  
  return breaks;
}
