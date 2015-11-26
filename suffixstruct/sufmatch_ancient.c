
/*
 *  sufmatch.c
 *  functions for matching in suffixarrays
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/19/06 15:13:18 CET
 *  
 */


 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "mathematics.h"
 #include "sufarray.h"
 #include "sufmatch.h"
 #include "mm.h"
 #include "intsequence.h"
 #include "list.h"
 #include "depictseqs.h"
 #include "gnuplot_i.h"
/*------------------------------- sufSubstring -------------------------------
 *    
 * retrieve the longest substring matches in a suffixarray
 * 
 */
 
PairSint*
sufSubstring (void *space, Suffixarray *arr, Uint *pattern, 
			  Uint len, Uint sublen)
{
  	Uint i;
	PairSint *res, d;

	if (len <= sublen) 
	{
	  return NULL;
	}
	
	res = ALLOCMEMORY(space, NULL, PairSint, len-sublen);
	for(i=0; i < len-sublen; i++) 
	{
		d=mmsearch(arr, &pattern[i], sublen, 0, 0, arr->numofsuffixes-1);
		res[i].a=d.a;
	  	res[i].b=d.b;	  
	}
	
	return res;
}


/*------------------------------ reportSufmatch ------------------------------
 *    
 * returns a beautified match string
 * 
 */
 
void
reportSufmatch (Suffixarray *a, PairSint *matches, 
				Uint len, Uint threshold,
				IntSequence **s)
{
    Uint i, j, idx;
	char *info;

	for (i=0; i < len; i++) {
	   if (matches[i].b >= ((matches[i].a)+threshold)) {
		   /*valid matches*/
		  for (j=matches[i].a; j <= matches[i].b; j++) {	
             idx = getMultiSeqIndex(a->seq, a->suffixptr[a->suftab[j]]);
			 info = s[idx]->description;
			 printf("[%d]:\t %s\n", j, info);
		  }
	   }
    }
	
    return ;
}


/*---------------------------------- cmp_* -----------------------------------
 *    
 * compare functions for clibs's qsort and bsearch
 * 1. cmp_suffixno : used in rankSufmatch to sort sequence numbers
 * 2. cmp_ranks    : used in rankSufmatch to sort sequence ranks
 * 3. cmp_ranks_ptr: used in rankSufmatchList to sort sequence ranks
 * 
 */

int
cmp_score(const void *a, const void *b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
    double frac_first, frac_second;

	if (first->score == 0 && second->score == 0) {
		if(first->count > second->count) return 1;
		if(first->count < second->count) return -1;
	}
	
	frac_first = (double) first->score ;//first->count;
	frac_second = (double) second->score ;//second->count;

	if(frac_first > frac_second) return 1;
	if(frac_first < frac_second) return -1;
	
/*	if(first->score > second->score) return 1;
	if(first->score < second->score) return -1;
*/
	return 0;
}


int 
cmp_suffixno(const void *a, const void* b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
 
	
	if(first->id > second->id) return 1;
	if(first->id < second->id) return -1;

	return 0;
}


int 
cmp_ranks(const void *a, const void* b) {
    Matchtype *first = (Matchtype*)a;
    Matchtype *second =(Matchtype*)b;	
  
	if(first->count > second->count) return 1;
	if(first->count < second->count) return -1;

	return 0;
}

int 
cmp_rank_ptr(const void *a, const void* b) {
    Matchtype **first =  (Matchtype**)a;
    Matchtype **second = (Matchtype**)b;	
  
	if(first[0]->count > second[0]->count) return 1;
	if(first[0]->count < second[0]->count) return -1;

	return 0;
}




/*------------------------------ freeMatchtype -------------------------------
 *    
 * a function to delete Matchtypes in a list
 * 
 */
 
void
freeMatchtype (void *space, void *data)
{
    Matchtype *d = (Matchtype*)data;

	FREEMEMORY(space, d->pos);
	FREEMEMORY(space, data);
}


/*------------------------------- rankSufmatchFct ------------------------------
 *    
 * ranks matches with a suffix array
 * given in an array of PairSint of length len. Sorting is done
 * by several calls to clib's qsort. For each item of the sorted
 * array a handler-function is invoked.
 * 
 */

void
rankSufmatchFct(void *space, Suffixarray *a, PairSint *matches, 
			 Uint len, Uint threshold, Uint K, Uint S, IntSequence **s,
			 void (*handler)(void *, Matchtype *, IntSequence **, int, void *), 
			 IntSequence *queryseq, void *info)
{
	Matchtype key, *cur, 
			  *occ=NULL;
	void *r;	
	Sint i, j, l, k=0;
	Uint *ptr;
	Uint temp;
	char *constr;
	Uint *consensus;
	double *consensus_double;
	IntSequence *consensusSequence;
	gnuplot_ctrl *h;
		
	consensus = ALLOCMEMORY(space, NULL, Uint, len);
	memset(consensus, 0, sizeof(Uint)*len);
	
	for(i=0; i < len; i++) {
		if(matches[i].b >= ((matches[i].a))) {
		
		  	for(j=matches[i].a; j <= matches[i].b; j++) {
		        key.id = getMultiSeqIndex(a->seq, a->suffixptr[a->suftab[j]]);	
				r=bsearch(&key, occ, k, sizeof(Matchtype), cmp_suffixno);
			
				if (r == NULL) {
				    occ = ALLOCMEMORY(space, occ, Matchtype, k+1);
				    occ[k].id = key.id;
					occ[k].count = 1;
					occ[k].pos = ALLOCMEMORY(space, NULL, Uint, 1);		
					occ[k].pos[0] = i;
					occ[k].score = 0;
					/*score the matches if no < K*/
					if ((matches[i].b-matches[i].a) < K) {
					    ptr = (a->suffixptr[a->suftab[j]]);
						for (l=0; l < S; l++){
							temp = (Uint) *ptr;
						    occ[k].score += (float)(temp/(1+matches[i].b-matches[i].a));
							ptr++;
						}
						consensus[i] ++;
					}			
		
					k++; 
					qsort(occ, k, sizeof(Matchtype), cmp_suffixno);
				} else {
					cur=((Matchtype*)r);
					cur->count++;
					cur->pos = ALLOCMEMORY(space, cur->pos, Uint, cur->count);
					cur->pos[(cur->count)-1]=i;
					/*score matches if no < K*/
					if ((matches[i].b-matches[i].a) < K) {
					  ptr = (a->suffixptr[a->suftab[j]]);
					  for (l=0; l < S; l++) {
						temp = (Uint) *ptr;
						cur->score += (float) (temp/(1+matches[i].b-matches[i].a));
						ptr++;
					  }	
					  consensus[i] ++;
					}
				}
			}
		}
	}
	
	qsort(occ, k, sizeof(Matchtype), cmp_score);
	l=0;
	for (i=k; i > 0 && l<=100; i--) {
	    if (occ[i-1].count > threshold) {
		    l++;
		    printf("IMBISS-Match %d: \n", l);
			handler(space, &occ[i-1], s, len, info);	
		}
	    FREEMEMORY(space, occ[i-1].pos);
	}
    

    for (i=0; i < len; i++) {
	    temp = 0;
		for(l=0; l < S; l++) {
			temp += queryseq->sequence[i+l]-1;
		}
		consensus[i] *= (temp/S);
	}
	
	consensusSequence = initSequence(space);
	consensusSequence->sequence = consensus;
	consensusSequence->length = len;
	constr = printSequence(space, consensusSequence,60);
	printf("%s\n", constr);

	consensus_double=ALLOCMEMORY(space, NULL, double, len);

	h = gnuplot_init();
	gnuplot_setstyle(h, "lines");
	for(i=0; i < len; i++) {
	  consensus_double[i] = (double) consensus[i];
	}
    gnuplot_plot_x(h, consensus_double, len, "blaPlot");
	
	FREEMEMORY(space, occ);    
	FREEMEMORY(space, constr);
	FREEMEMORY(space, consensus);
	FREEMEMORY(space, consensus_double);
	FREEMEMORY(space, consensusSequence);
}

/*------------------------------- rankSufmatch -------------------------------
 *    
 * ranks matches with a suffix array
 * given in an array of PairSint of length len. Sorting is done
 * by several calls to clib's qsort.
 * 
 */

void
rankSufmatch(void *space, Suffixarray *a, PairSint *matches, 
			 Uint len, Uint threshold, IntSequence **s)
{
	Matchtype key, *cur, 
			  *occ=NULL;
	void *r;
	char *seq;
	Sint i, j, k=0;
		
	for(i=0; i < len; i++) {
		if(matches[i].b >= ((matches[i].a))) {
		
		  	for(j=matches[i].a; j <= matches[i].b; j++) {
		        key.id = getMultiSeqIndex(a->seq, a->suffixptr[a->suftab[j]]);	
				r=bsearch(&key, occ, k, sizeof(Matchtype), cmp_suffixno);
			
				if (r == NULL) {
				    occ = ALLOCMEMORY(space, occ, Matchtype, k+1);
				    occ[k].id = key.id;
					occ[k].count = 1;
					occ[k].pos = ALLOCMEMORY(space, NULL, Uint, 1);		
					occ[k].pos[0] = i;
					occ[k].score = (double) s[key.id]->sequence[i]/len;
					k++; 
					qsort(occ, k, sizeof(Matchtype), cmp_suffixno);
				} else {
					cur=((Matchtype*)r);
					cur->count++;
					cur->pos = ALLOCMEMORY(space, cur->pos, Uint, cur->count);
					cur->pos[(cur->count)-1]=i;
					cur->score = (double) s[cur->id]->sequence[i]/len;
				}
			}
		}
	}
	
	qsort(occ, k, sizeof(Matchtype), cmp_ranks);
	
	for (i=k; i > 0; i--) {
	    if (occ[i-1].count > threshold) {
	    printf("%d\t%s\t%d\t", occ[i-1].id, s[occ[i-1].id]->url, 
		 				  occ[i-1].count);
	    seq= depictSequence(space, len, 20, occ[i-1].pos, occ[i-1].count,'*');
		printf("[%s]\n", seq);
		printf("%s\n", s[occ[i-1].id]->description);
		
		FREEMEMORY(space, seq);
		}
	    FREEMEMORY(space, occ[i-1].pos);
	}
	FREEMEMORY(space, occ);
	
	printf("found %d (partial) matches of %d.\n", k, len);
    
}



