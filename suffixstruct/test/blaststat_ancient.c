
/*
 *  blaststat.c
 *  an implementation of blast statistics
 *  completely re-written: no fancy shit implementation
 *  
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 02/26/07 21:18:34 CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mathematics.h"
#include "memory.h"
#include "intsequence.h"
#include "falphabet.h"
#include "basic-types.h"
#include "blaststat.h"

#define BLAST_SUM_LIMIT 0.001 /*speedy*/
#define BLAST_ITER_MAX  100

/*--------------------------------  seqfreq  ---------------------------------
 *    
 * returns the frequency of query characters in alphabet A
 * Function accepts a query sequence and the corresponding
 * alphabet. Frequencies are returned in a vector.
 * 
 */
 
double *
seqfreq (void *space, IntSequence *seq, FAlphabet *alpha)
{
    Uint l, i;
	double *f, sum =0; 
	  
	l = alpha->domainsize;
	f = ALLOCMEMORY(space, NULL, double, l);
	memset(f, 0, sizeof(double)*l);
	for(i=0; i < seq->length; i++) {
		f[seq->sequence[i]]++;
		sum++;
	}

	for(i=0; i < l; i++) {
		f[i]/=sum;
	}
	
	return f;
}


/*---------------------------------- dbfreq ----------------------------------
 *    
 * returns the frequencies in a set of sequences for each char in A 
 * ResFreqNormalize
 */
 
double*
dbfreq (void *space, IntSequence **set, Uint setsize, FAlphabet *alpha,
		double norm)
{
  	Uint i, j, l, sum=0;
	double *f;
	
	l = alpha->domainsize;
	f = ALLOCMEMORY(space, NULL, double, l);
	memset(f, 0, sizeof(double)*l);
	
	for(i=0; i < setsize; i++) {
		for(j=0; j < set[i]->length; j++) {
			f[set[i]->sequence[j]]++;
			sum++;
		}
	}

	for(i=0; i < l; i++) {
		f[i] /= sum;
		f[i] *= norm;
	}

	
	return f;
}


#define BASIC_PROBABILITY 0.000001

/*-------------------------------- logoddscr ---------------------------------
 *    
 * calculations of log-odds scores (Karlin Altschul PNAS 89, example iii)
 * 
 */

double*
logoddscr (void* space, double* dbf, double *sf, FAlphabet *alpha)
{
	 Uint i;
	 double temp;
	 double *logodd=NULL;
 
	 logodd = ALLOCMEMORY(space, NULL, double, alpha->domainsize);
	 memset(logodd, 0, sizeof(double)*alpha->domainsize);

	 for(i=0; i < alpha->domainsize; i++) {
		if (dbf[i] > 0) {
		    temp = (sf[i]+0.0000000001)/dbf[i];
			logodd[i] = log(temp)*10;	
		} else {
		  	logodd[i] = (double)0;
		}
	 }
	 
	return logodd;
}


void writescores(double *scr, Uint l, char* fname) {
	FILE *fp;
	Uint i;
	
	fp = fopen(fname, "w");
	if (fp == NULL) {
	    fprintf(stderr, "error: unable to open file.\n");
		return;
	}
	
	for(i=0; i < l; i++) {
		fprintf(fp, "%f\n", scr[i]);
	}

	fclose(fp);
}

/*-------------------------------- scorefreq ---------------------------------
 *    
 * calculates the frequencies of all given integer scores. 
 * Function accepts an array of scores and arrays containing
 * query and database frequencies for each char in alphabet A.
 * returns an array of frequencies. Additionally avg is filled.
 * 
 */
 
double *
scorefreq ( void *space, double *scr, Uint noofscr, FAlphabet *alpha, 
			double *sf, double *dbf)
{
	Uint i;
	double sum = 0;	
	double *f;

	f = ALLOCMEMORY(space, NULL, double, noofscr);
	memset(f, 0, sizeof(double)*noofscr);
	
	for(i=0; i < noofscr; i++) {
		f[i] += sf[i]*dbf[i];
		sum += sf[i]*dbf[i];
	}
	
	
	if (sum > 0.0001) {
		for(i=0; i < noofscr; i++) {
			f[i] /= sum;
		}
	}

  	return f;
}




double
checklambda (double *scr, Uint noofscores, double* sf, double avg, double e_lamb){
	Uint i;
    double val=0;
	
  	for(i=0; i < noofscores; i++) {
		val += sf[i]*exp((e_lamb*scr[i])); 
	}

	return val;
} 
		


/*-------------------------------- score_evd ---------------------------------
 *    
 * extreme value distribution for given scores (supplied in info)
 * 
 */
 
double
score_evd (double lambda, void* info)
{
  	Uint i;
	double sum=0;
	evdparam *p =(evdparam*) info;
	
	for(i=0; i < p->noofscores; i++) {
		sum += p->probs[i]*exp((double) lambda*p->scores[i]);
	}
	
	return (sum-1.0);
}



/*---------------------------------- lambda ----------------------------------
 *    
 * returns the Karlin-Altschul lambda given scores and frequencies
 * 
 */
 
double
lambda (Uint *sort, double *scr, Uint noofscores, double* sf, 
		double avg, double initlambda, double tol)
{
   	Uint i, j, step=1,
		 iterations=100,
		 maxNewton=100+17,
		 wasNewton=0,
		 isNewton=0;
	double e_init,
		   e_lambda,
		   f=4, f_prime,
		   p, y,
		   a=0,
		   b=1,
		   fold;
	
    if (avg >= 0) {
	  return -1.0;
	}

	/*implement the gcd method to improve integration here*/
    /*stepsize d set to 1, that is the smallest cd for int*/
	e_init = exp(-initlambda);
	e_lambda = ( 0 < e_init && e_init < 1) ? e_init : 0.5;


	for (i=0; i < iterations; i++) {
	  	wasNewton = isNewton;
		isNewton = 0;
		printf("iter: %d, lambda: %19.16e, %19.16e \n", i, e_lambda, e_init);	
		f = scr[sort[0]];
		f_prime = 0;
		
		/*iter the negative scores*/
		for(j=step; j < noofscores && scr[sort[j]] < 0; j+=step) {
			f_prime = e_lambda * f_prime + f;
			f = e_lambda * f + scr[sort[j]];
		}
		
		/*should always be true*/
		/*while (scr[sort[j]] == 0) {*/
			f_prime = e_lambda * f_prime +f;
			f = e_lambda * f + scr[sort[j]] -1;
		/*	j += step;
		} */

		for( ; j < noofscores && scr[sort[j]] >= 0; j+= step) {
			f_prime = e_lambda * f_prime + f;
			f = e_lambda * f + scr[sort[j]];
		}

		printf("lambda: %19.16e f:%f g:%f \n", e_lambda, f, f_prime);
		if (f > 0) {
			a = e_lambda;
		} else if (f < 0) {
			b = e_lambda;
		} else {
			break;
		}

		if(b - a < 2*a * (1-b)*tol) {
			e_lambda = (a+b) / 2; break;
		}

		if (i >= maxNewton || f_prime >= 0 ||
			((wasNewton) && fabs(f) > 0.9 * fabs(fold)) )
		{
			e_lambda = (a+b)/2;
		} else {
			p = - f/f_prime;
			y = e_lambda + p;
			if (y <= a || y >= b) {
				e_lambda = (a+b)/2;
			} else {
				isNewton = 1;
				e_lambda = y;
				if (fabs(p) < e_lambda * (1-e_lambda) * tol)
				  break;
			}
		}
	}

	printf("e: %19.16e\n", e_lambda);
	return (-log(e_lambda)/step);
}


/*-------------------------------- relentropy --------------------------------
 *    
 * calculates the relative entropy of probabilities and target
 * frequencies.
 *
 * basics: 
 * scr is sorted - the observed minimum score is in scr[0]
 * probs contains the probabilities of scrs - eg. q[scr[0]]
 * 
 */
 
double
relentropy (Uint *sort, double *scr, Uint noofscr, double *q, double lambda)
{
  	double e_lambda,
		   sum=0, 
		   scale, H, av,
		   sumval, K, Ktmp;
	double *p, *P, *ptrP, *ptr1, *ptr2, *ptr1e, *resc, *scr2;
	int i, j, first, last, hi, lo, range;
	
	
  	if (lambda < 0)
	  return -1;

    /*sum the probabilities*/	
	for (i=0; i < noofscr; i++) {
		sum += q[i];
	}	

	/*normalize probabilites*/
	p = ALLOCMEMORY(space, NULL, double, noofscr);
	for (i=0; i < noofscr; i++) {
		p[i] = q[i] / sum;
	}

	/*resort*/
	resc  = ALLOCMEMORY(space, NULL, double, noofscr);
	scr2  = ALLOCMEMORY(space, NULL, double, noofscr);
	for(i=0; i < noofscr; i++) {
		resc[i] = p[sort[i]];
		scr2[i] = scr[sort[i]];
	}
	FREEMEMORY(space, p);
	p = resc;

	/*get average*/
	for (av=0.0, i=0; i < noofscr; i++) {
		av += p[i]*(scr2[i]+1)*(lambda*(scr2[i]+1));
	}

	/*get P*/
	P = ALLOCMEMORY(space, NULL, double, (150 * noofscr +1));
	
	sumval=0;
	lo=0;
	hi=0;
	
	/*fancy sum /frac{1}{k=j} E[e^{\lambda*S_k}; S_k <0]+Prob(S_k>=0) }*/
	/*			^frac here		^calculated at the end of the loop     */
	
	for (*P = 1.0, sum = 1.0, j = 1;
			j <= 150 && sum > 0.00001; sumval += sum /= j++) {
	
	  first=last=noofscr-1;

	  for (ptrP = P + (hi += scr2[noofscr-1]) - (lo += scr2[0]);
		  ptrP >= P; *ptrP-- = sum) {
	  	ptr1   = ptrP - first;
		ptr1e  = ptrP - last;
		ptr2   = p + first;

		for(sum = 0.0; ptr1 >= ptr1e;) {
			sum += *ptr1-- * *ptr2++;
		}
		if (first) {
			--first;
		}
		if (((int) (ptrP-P)) <= range) {
			--last;
		}
	  }
	  
	  /*get \sum E[e^{\lambda*S_{k}}, S_k < 0]*/
	  for(sum=0.0, i = lo; scr2[i] < 0; ++i) {
	  	sum += *++ptrP * exp(lambda*scr2[i]);
	  }
	  for(; i <= hi; ++i) {
	  	sum += *++ptrP;
	  } 
	}
	Ktmp = (double) (exp(-2 * sumval));
	K = Ktmp / (av * (1.0 - exp(-lambda)));

	/*
	e_lambda = exp(-lambda);
	sum = scr[sort[0]] * q[sort[0]];

	for(i=1; i < noofscr; i++){
		sum = scr[sort[i]] * q[sort[i]] + e_lambda * sum;
	}

	scale = power(e_lambda, scr[sort[i-1]]);
	if (scale > 0) {
		H = (lambda * sum/scale);
	} else {
		H = (lambda * exp (lambda * scr[sort[i-1]] + log(sum)));
	}
	
	return H;*/
	return K;
}



double
C_star (Uint *sort, double *scr, Uint noofscr, double avg, double *q, 
		  double lambda, double H) {

return 0;
}


/*--------------------------------- entropyK ---------------------------------
 *    
 * calculate Karlin-Altschul K (PNAS 87 (1990))
 * greatest common divisor set to 1. Fix?
 * 
 */
 
double
entropyK (Uint *sort, double *scr, Uint noofscr, double avg, double *q, 
		  double lambda, double H)
{

  	int i,j;
  	int first, last, range;
	int itermax;
	int lo, hi; /*alignment scores*/
	
	double *align;
	double *ptrP, *ptr1, *ptr2, *ptr1e;
	double sumlimit;
	double outer; /*outer sum*/
	double inner, oldinner, ancientinner; /*inner sums*/
	
	double K;
	double HdivLambda;
	double expminLambda;
	
    if (lambda <= 0 || H <= 0 || avg >=0) {
	  	fprintf(stderr, "Error: H and lambda expected to be < 0\n");
		return -1;
	}
	
	range = scr[sort[noofscr-1]]-scr[sort[0]];
	/*gcd?*/
	
	HdivLambda = H/lambda;
	expminLambda = exp((double) -lambda);
	
	/*first case*/
	if(scr[sort[0]]==-1 && scr[sort[noofscr-1]]==1){
		K = q[sort[0]] - q[sort[noofscr-1]];
	  	K *= K;
	  	K /= q[sort[0]]; 
	 	return K;
	}

	/*second case*/
	if (scr[sort[0]] == -1 || scr[sort[noofscr-1]] == 1) {
		if (scr[sort[noofscr-1]] != 1) {
			HdivLambda = (avg*avg)/HdivLambda;
		}
		return HdivLambda *(1.0 - expminLambda);
	}

	/*third case*/
	sumlimit = BLAST_SUM_LIMIT;
	itermax = BLAST_ITER_MAX;

	align = (double*) calloc ((itermax*range+1), sizeof(*align));
	if (align == NULL) {
		fprintf(stderr, "Error: Blast alloc failed.\n");
		exit(-1);
	}

	outer = 1;
	lo = hi = 0;
	align[0] = inner = oldinner = ancientinner = 1;
	
	for(i=0; i < itermax && inner > sumlimit; outer += inner /= ++i) {
		
		first = last = range;
		lo += scr[sort[0]];
		hi += scr[sort[noofscr-1]];

		for(ptrP = align + (hi - lo); ptrP >= align; *ptrP-- = inner) {	
		  	ptr1  = ptrP - first;
			ptr1e = ptrP - last;
			ptr2 = &q[sort[0]] + first; /*tricky?*/
			for(inner = 0; ptr1 >= ptr1e; ) {
			  inner += *ptr1-- * *ptr2++;
			}
			if (first) --first;
			if (ptrP - align <= range) --last;			
		}

		/*Hoernchen*/
		inner = *++ptrP;
		for(j = lo+1; j < 0; j++) {
			inner = *++ptrP + inner * expminLambda;
		}
		inner *= expminLambda;

		for(;j <= hi; ++j) {
			inner += *++ptrP;
		}

		ancientinner = oldinner;
		oldinner = inner;
	}
	printf("\n");
	K = -exp((double)-2.0*outer) / (HdivLambda*BLAST_Expm1(-(double)lambda));
	free(align);

	return K;
}

