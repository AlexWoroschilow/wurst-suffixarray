
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
  
	if(first->score > second->score) return 1;
	if(first->score < second->score) return -1;

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



/*----------------------------- rankSufmatchList -----------------------------
 *    
 * ranks matches with a suffix array
 * given in an array of PairSint of length len. Sorting is done
 * using a linearly linked list.
 * 
 */

void
rankSufmatchList (void *space, Suffixarray *a, PairSint *matches, 
			 	  Uint len, Uint threshold, IntSequence **s)
{
	List *l=NULL; 
	Matchtype *occ=NULL, **m;
	
	char *seq;
	Sint i, j;
	Uint id, cur;
		
    for (i=0; i < len; i++) {
	   if (matches[i].b >= ((matches[i].a)+threshold)) {
		  for (j=matches[i].a; j <= matches[i].b; j++) {
		     id=getMultiSeqIndex(a->seq, a->suffixptr[a->suftab[j]]);	
			 if (l != NULL) {
			    cur=l->firstNode; 
				while (((Matchtype*)l->nodes[cur].data)->id < id && 
					  				l->nodes[cur].next != LISTNILVALUE) {
				   cur=l->nodes[cur].next;
				}
			 } 
				
			 if (l == NULL ||((Matchtype*)l->nodes[cur].data)->id != id) {
			    occ = ALLOCMEMORY(space, NULL, Matchtype, 1);
				occ->id = id;
				occ->count = 1;
				occ->pos = ALLOCMEMORY(space, NULL, Uint, 1);		
				occ->pos[0] = i;
					
				if (l == NULL) {
				   l=initList(space, 10);
				   insertBefore(space, l, LISTNILVALUE, occ);
				  
				} else {
					
				   if (((Matchtype*)l->nodes[cur].data)->id > id) 
					  insertBefore(space, l, cur, occ);
				   else 
					  insertAfter(space, l, cur, occ);		
				}
			 } else {
		        occ = ((Matchtype*)l->nodes[cur].data);
				occ->count++;
				occ->pos = ALLOCMEMORY(space, occ->pos, Uint, occ->count);	
				occ->pos[occ->count-1] = i; 
			 }
	      }
	   }
	}

	m = (Matchtype**) dataList(space, l);
	qsort(m, l->length, sizeof(Matchtype*), cmp_rank_ptr);
	
	for (i=l->length; i > 0; i--) {
		printf("%d\t%s\t%d", m[i-1]->id, s[m[i-1]->id]->description,
			                 m[i-1]->count);
		seq = depictSequence(space, len, 20, m[i-1]->pos, m[i-1]->count, '*');
		printf("\t[%s]\n", seq);
		FREEMEMORY(space, seq);
		FREEMEMORY(space, m[i-1]->pos);
		FREEMEMORY(space, m[i-1]);
	}
		
	printf("found %d (partial) matches of %d.\n", l->length, len);
	FREEMEMORY(space, m);
	wrapList(space, l, NULL);	
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
					k++; 
					qsort(occ, k, sizeof(Matchtype), cmp_suffixno);
				} else {
					cur=((Matchtype*)r);
					cur->count++;
					cur->pos = ALLOCMEMORY(space, cur->pos, Uint, cur->count);
					cur->pos[(cur->count)-1]=i;
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
			 Uint len, Uint threshold, IntSequence **s,
			 void (*handler)(void *, Matchtype *, IntSequence **, int, void *), 
			 void *info)
{
	Matchtype key, *cur, 
			  *occ=NULL;
	void *r;	
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
					occ[k].score = (float)s[key.id]->sequence[i]/len;
					k++; 
					qsort(occ, k, sizeof(Matchtype), cmp_suffixno);
				} else {
					cur=((Matchtype*)r);
					cur->count++;
					cur->pos = ALLOCMEMORY(space, cur->pos, Uint, cur->count);
					cur->pos[(cur->count)-1]=i;
					cur->score += (float) s[cur->id]->sequence[i] /len;
				}
			}
		}
	}
	
	qsort(occ, k, sizeof(Matchtype), cmp_score);
	
	for (i=k; i > 0; i--) {
	    if (occ[i-1].count > threshold) {
			handler(space, &occ[i-1], s, len, info);	
		}
	    FREEMEMORY(space, occ[i-1].pos);
	}
	FREEMEMORY(space, occ);    
}


