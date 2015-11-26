
/*
 * mathematics.c
 * implemtation of various mathematical functions
 *
 * @author Steve Hoffmann
 * @date Wed 22 Nov 2006
 *
 */

#include "mathematics.h"
#include <math.h>

 void *initArray(void *space, int size, Uint datatype) {
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


void dumpMatrix_int(int *M, int m, int n) {
	int i,j;

	for (i=0; i < m; i++) {
		for (j=0; j < n; j++){
			printf("%d ", MATRIX2D(M,n,i,j));
		}
			printf("\n");
	}
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

	for (j=i+1; j < v->length; j++ )
		if(e[i-1]>=e[j]) break;
	
	SWAPVEC(i-1, j-1, v);
	REVERSEVEC(i, (v->length)-1, v);

	return 1;
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
  	double y;

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

	y = 1.;
	while(n > 0) {
		if (n & 1) {
			y *= x;
		}
		n /= 2;
		x *= x;
	}

	return y;
}


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
uniroot (double start, double end, double (*f)(double, void*), double tolx, void* info)
{	
  	double a, b, c;
	double fa, fb, fc;
	double prev;
	double currenttol;
	double p, q, new_step;
	double cb, t1, t2;

	a = start; b= end; fa = (*f)(a,info); fb=(*f)(b,info);
	c = a; fc = fa;
	
	if ((fa > (double) 0 && fb > (double) 0) || (fa < (double)0 && fb < (double)0)) {
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


double*
transpose (double *a, Uint m, Uint n) {
	double *t;
	double 	r;
	int		i,
			j=0,
			k=0;

	t = (double*) INITMATRIX2D(t, n, m, sizeof(double));

	for(i=0; i < m*n; i++) {
	  if(i % n) { j++; k=0;} 	
	  MATRIX2D(t, n, k, j) = a[i];
	  k++;
	}

	return t;
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


