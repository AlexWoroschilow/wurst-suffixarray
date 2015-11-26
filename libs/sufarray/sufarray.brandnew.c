
/*
 *  sufarray.c
 *  implementations for enhanced suffix arrays
 *  for large integer alphabets
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/11/06 14:56:57 CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "sufarray.h"
#include "charsequence.h"
#include "falphabet.h"
#include "stack.h"
#include "sort.h"
#include "list.h"
#include "vtprogressbar.h"
#include "aluruSort.h"
#include "queuedef.h"
#include "debug.h"
#include "bitArray.h"


void
destructinterval(void *space, void *data) {
  FREEMEMORY(space, (PairUint*) data);
}

/* ------------------------------ cmpCharSequence ----------------------------
 *    
 * function to compare CharSequences for mulitkey sort (sort.c)
 * 
 */

  Uint
cmpCharSequence (Uint a, Uint b, Uint depth, void *data, void *info)
{
  char **s = (char**) data;	
  Uint *end;

  /*quick fix to meet end of multiintsequence criterion*/
  if (info == NULL) {
    if(s[b][0] == (char) 127) {
      if (s[a][depth] == (char) 127) {
        return 0;
      }
      return 1;
    }
  } else {
    end = (Uint*) info;
    if (*end == b) {
      if (s[a][depth] == (char) 127) {
        return 0;
      }
      return 1;
    }
  }


  /*real comparison*/
  if (s[a][depth] > s[b][depth]) return 1;
  if (s[a][depth] < s[b][depth]) return 2;

  return 0;
}


/*------------------------------- getSuffixPtr -------------------------------
 *    
 * returns an array of pointers for all "real" suffixes in mseq.
 * 
 */

  char**
getSuffixPtr(void *space, MultiCharSeq *mseq, char delim, char sentinel) 
{
  Uint i,r,k=0;
  char **ptr;

  ptr = ALLOCMEMORY(space, NULL, char*, 
      (mseq->totallength)-(mseq->numofsequences)+1);

  for(i=0; i < mseq->totallength; i++) {
    if((r=mseq->sequences[i]) != delim){
      ptr[k]= &mseq->sequences[i];	  
      k++;
    }
  }

  return ptr;
}


/* ---------------------------- constructSufArr -----------------------------
 *    
 * constructs a suffix array from an (unsigned) integer sequence
 * should be working in O(n log n). It uses a multi-key quicksort.
 * 
 */

  Suffixarray*
constructSufArr(void *space, 
    CharSequence **s, 
    Uint len, 
    FAlphabet* alphabet)
{

  Uint i, numofsuffixes,
  *sorted, 
  *inv_suftab;
  char **suffixptr;
  MultiCharSeq *mseq;
  char save;
  Suffixarray *arr;

  mseq = concatCharSequences(space, s, len, (char)254, (char)127);
  numofsuffixes = (mseq->totallength - mseq->numofsequences)+1; 

  fprintf(stderr, "alphabet of size (%d): %s\n", mseq->mapsize, mseq->map);
  fprintf(stderr, "allocating space for %d suffixes (%d bytes)\n", numofsuffixes, 2*numofsuffixes*sizeof(Uint));
  inv_suftab = ALLOCMEMORY(space, NULL, Uint , numofsuffixes);
  arr = ALLOCMEMORY(space, NULL, Suffixarray, 1);

  fprintf(stderr, "constructing suftab.\n");
  suffixptr = getSuffixPtr(space, mseq, (char)254, (char)127);
  save = mseq->sequences[numofsuffixes-1]; 

  //mseq->sequences[numofsuffixes-1]=0; 
  sorted = alurusort(space, mseq->sequences, &(numofsuffixes));
  //  mseq->sequences[numofsuffixes-1]=save;

  /*  sorted = quickSortMultikey (space, suffixptr, numofsuffixes, 
      cmpCharSequence, numofsuffixes-1, NULL);     
      */
  fprintf(stderr, "constructing inv_suftab.\n");
  for (i=0; i < numofsuffixes; i++) {
    inv_suftab[sorted[i]]=i;
  }

  arr->seq = mseq;
  arr->numofsuffixes = numofsuffixes;
  arr->suffixptr = suffixptr;
  arr->suftab = sorted;
  arr->inv_suftab = inv_suftab;

  return arr;
}


void
writeSuffixarray(Suffixarray *s, char *filename) {
  FILE *fp; 
  size_t nmemb;

  fp = fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr,"Couldn't open file. Exit forced.\n");
    exit(-1);
  }

  nmemb = (size_t) s->numofsuffixes;
  fwrite(&nmemb, sizeof(size_t), 1, fp);
  fwrite(s->suftab, sizeof(Uint), nmemb, fp);

  fclose(fp);

}


void
readSuffixarray(Suffixarray *s, char *filename) {
  FILE *fp; 
  size_t nmemb;
  Uint i;
  Uint *suftab;

  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr,"Couldn't open file. Exit forced.\n");
    exit(-1);
  }


  fread(&nmemb, sizeof(size_t), 1, fp);
  suftab = ALLOCMEMORY(NULL, NULL, Uint, nmemb);
  fread(suftab, sizeof(Uint), nmemb, fp);

  for(i=0; i < nmemb; i++) {
    if(suftab[i] != s->suftab[i]) exit(-1);
  }
  printf("suftab successfully read nmemb %d, num %d.\n", nmemb, s->numofsuffixes);

  fclose(fp);

}



/*------------------------------ destructSufArr ------------------------------
 *    
 * destruct a suffix array.
 * 
 */
  void
destructSufArr (void *space, Suffixarray *arr)
{
  FREEMEMORY(space, arr->suftab);
  if (arr->lcpctab != NULL)
    FREEMEMORY(space, arr->lcpctab);
  if (arr->inv_suftab != NULL)
    FREEMEMORY(space, arr->inv_suftab);
  if (arr->suffixptr != NULL)
    FREEMEMORY(space, arr->suffixptr);
  if (arr->seq != NULL)
    destructMultiCharSeq(space, arr->seq);
  FREEMEMORY(space, arr);

  return ;
}


inline Uint
lcp(Suffixarray *s, Uint i) {
  PairUint *ret;
  Uint val;

  /*  return s->lcptab[i];*/

  if(s->lcpctab[i] < 254) {
    return (Uint) s->lcpctab[i];
  } else { 
    val = i;
    ret=bsearch(&val, s->llvtab, s->llvcnt, sizeof(PairUint), cmp_PairUint_bsearch);
  }
  if (ret == NULL) {
    printf("lcp not found. Exit forced.\n");
    exit(-1);
  }

  return ret->b;
}


inline Uint
maxlcp(Suffixarray *s) {
  Uint i;
  Uint max = 0;

  for(i=0; i < s->numofsuffixes; i++) {
    if (lcp(s,i) > max) max = lcp(s,i);
  }
  return max;
}

/*------------------------------ computeLcpTab -------------------------------
 *    
 * computes the lcp tab from suftab and inv_suftab in O(n).
 * 
 */

  void
constructLcp (void *space, Suffixarray *arr)
{
  Uint i, j, k;
  Uint s, t, max = 0;
  int l=0;

  /*arr->lcptab = ALLOCMEMORY(space, NULL, Uint, arr->numofsuffixes);
    memset(arr->lcptab, 0, sizeof(Uint)*arr->numofsuffixes);*/

  arr->lcpctab = ALLOCMEMORY(space, NULL, unsigned char, arr->numofsuffixes);
  arr->llvcnt = 0;
  arr->llvtab = NULL;

  for(i=0; i < arr->numofsuffixes; i++) {    
    j = arr->inv_suftab[i];

    if (j > 0) {
      k = arr->suftab[j-1];
      s = arr->suffixptr[k]-arr->seq->sequences;
      t = arr->suffixptr[i]-arr->seq->sequences;


      /* j=i+1;
         s = arr->suftab[i];
         t = arr->suftab[j];
         */
      l=l-1;
      if (l < 0) l=0;

      while (arr->seq->sequences[s+l] == arr->seq->sequences[t+l]){ 
        l++;
      }

      /*    arr->lcptab[j] = l;*/
      if (l > max) max = l;
      if (l < 254) {
        arr->lcpctab[j] = (char) l;
      } else {
        arr->lcpctab[j] = 254;
        arr->llvtab = ALLOCMEMORY(space, arr->llvtab, PairUint, arr->llvcnt+1);
        arr->llvtab[arr->llvcnt].a = j;
        arr->llvtab[arr->llvcnt].b = l;
        arr->llvcnt++;
      }
    }
  }

  qsort(arr->llvtab, arr->llvcnt, sizeof(PairUint), cmp_PairUint_qsort);
  arr->maxlcp=max;
  arr->lcpctab[0]=0;

  return;
}

/*-------------------------------- dumpSufArr --------------------------------
 *    
 * dumps a suffix array to a screen
 * 
 */

  void
dumpSufArr (Suffixarray *arr)
{
  Uint i;

  for(i=0; i < arr->numofsuffixes; i++) {
    printf("%d \t %d \t %d \t %d \t %d \t %d\n", i, 
        arr->suftab[i], 
        lcp(arr,i),
        arr->inv_suftab[i], 
        *arr->suffixptr[arr->suftab[i]],
        arr->seq->sequences[arr->suffixptr[arr->suftab[i]]-arr->seq->sequences]);
  }

  return;
}

void
dumplcps(Suffixarray *arr) {
  Uint i, j, s, t;


  for(i=0; i < arr->numofsuffixes; i++) {
    if (lcp(arr,i) > 0) {
      s = arr->suffixptr[arr->suftab[i-1]]-arr->seq->sequences;
      t = arr->suffixptr[arr->suftab[i]]-arr->seq->sequences;
      printf("lcp of suffix %d and %d has length %d\t:\n", i-1, i, lcp(arr,i));
      for(j=0; j <= lcp(arr,i); j++) printf(" %d ", arr->seq->sequences[s+j]);
      printf("\n");
      for(j=0; j <= lcp(arr,i); j++) printf(" %d ", arr->seq->sequences[t+j]);
      printf("\n");
    }
  }
}

inline int
id(Suffixarray *s, Uint i) {
  PairSint *ret;
  int val;
  
  if(s->id[i] != (signed char)-128) {
    return (int) s->id[i];
  } else { 
    val = (int) i;
    ret = bsearch(&val, s->idvtab, s->idvcnt, sizeof(PairSint), cmp_PairSint_bsearch);
  }
  if (ret == NULL) {
    printf("id not found. Exit forced.\n");
    exit(-1);
  }
  return ret->b;
}


inline unsigned char
isnextlIndex(Suffixarray *s, Uint i) {
  return (lcp(s,s->chldtab[i].val) == lcp(s,i));
}

inline unsigned char
isdownIndex(Suffixarray *s, Uint i) {
  return (lcp(s,s->chldtab[i].val) > lcp(s,i));
}

inline unsigned char
isupIndex(Suffixarray *s, Uint i) {
  return (lcp(s,i) > lcp(s, i+1));
}

inline unsigned char
isLeaf(Suffixarray *s, Uint i, Uint j) {

  Uint     i1,
           i2;
  unsigned char notRoot;

  if (i==j) return 1;
  notRoot = (i > 0 || j < s->numofsuffixes-1);

  if(notRoot) {
    if (i < s->chldtab[j].val && s->chldtab[j].val <=j) {
      i1 = s->chldtab[j].val;
    } else {
      i1 = s->chldtab[i].val;
    }
    if (i !=  i1-1) return 0;
  } else {
    i1 = i;
  }

  while(isnextlIndex(s,i1) && !isupIndex(s,i1) 
      && s->chldtab[i1].val != 0) {
    i2 = s->chldtab[i1].val;
    if (i1 != i2-1) return 0;
    i1 = i2;
  }

  if(notRoot) {
    if (i1 != j) return 0;
  }

  return 1;
}
  

void
constructchildtab(void *space, Suffixarray *s) {
  Uint i;
  int lastIndex = -1;
  Stack *stack;

  s->chldtab = ALLOCMEMORY(space, NULL, childtab, s->numofsuffixes+1);
  memset(s->chldtab, 0, s->numofsuffixes*sizeof(childtab));
  stack = ALLOCMEMORY(space, NULL, Stack, 1);
  initStack(space, stack, 100000);

  stackpush(space, stack, 0);

  for(i=0; i < s->numofsuffixes; i++) 
  {
    while(lcp(s,i) < lcp(s,stacktop(stack))) {
      lastIndex = stackpop(stack);
      if(lcp(s,i) <= lcp(s,stacktop(stack)) && 
          lcp(s,stacktop(stack)) != lcp(s,lastIndex))
      {
        s->chldtab[stacktop(stack)].val  = lastIndex;
      }
    }
    if (lastIndex != -1) {
      s->chldtab[i-1].val = lastIndex;
      lastIndex = -1;
    }
    stackpush(space, stack, i);
  }

  /*construction of nextlIndex value*/
  destructStack(space, stack);
  initStack(space, stack, 10000);
  stackpush(space, stack,0);

  for(i=1; i < s->numofsuffixes; i++) {
    while(lcp(s,i) < lcp(s,stacktop(stack))) {
      stackpop(stack);
    }
    if (lcp(s,i) == lcp(s, stacktop(stack))) {
      lastIndex = stackpop(stack);
      s->chldtab[lastIndex].val = i;
    }
    stackpush(space, stack, i);
  }

  return;
}



/*------------------------------- computeID ----------------------------------
 *    
 * @brief performs a top down traversal on the tree represented 
 * by the suffix array to compute unique ids for each interval
 * @author Steve Hoffmann 
 *   
 */

void
computeId (void *space, Suffixarray *s) {
  Uint i, 
       a, 
       b, 
       l, 
       r;
  List* list;
  Queue q;


  emptyqueue(space, &q, 1000);  
  s->id = ALLOCMEMORY(space, NULL, char, s->numofsuffixes+2);
  memset(s->id, 0, sizeof(char)*s->numofsuffixes+2);

  s->idvtab = ALLOCMEMORY(space, NULL, PairSint, 1);
  s->idvcnt = 1;

  a = 0;
  b = s->numofsuffixes-1;

  s->id[0] = (signed char) -128;
  s->idvtab[0].a = 0;
  s->idvtab[0].b = (s->numofsuffixes-1);


  enqueue(space, &q, a);
  enqueue(space, &q, b);

  while (!queueisempty(&q)) {
        
    a = dequeue(&q);
    b = dequeue(&q);

    list = getChildintervals(space, s, a, b);  

    for (i=0; i < list->length; i++) {

      l = ((PairUint*)list->nodes[i].data)->a;
      r = ((PairUint*)list->nodes[i].data)->b;         
        
      if (l < r) {
        if(s->id[l] == 0) {
          if (r-l <= 127){
            s->id[l] = (signed char) r-l;
          } else {
            s->id[l] = (signed char)-128;
            s->idvtab = ALLOCMEMORY(space, s->idvtab, PairSint, s->idvcnt+1);
            s->idvtab[s->idvcnt].a = l;
            s->idvtab[s->idvcnt].b = r-l;
            s->idvcnt = s->idvcnt+1;
            // qsort(s->idvtab, s->idvcnt, sizeof(PairSint), cmp_PairSint_qsort);
          }

          /*if (id(s,l) != r-l) {
            printf("[%d,%d]=[%d,%d]", l, r, l, l+id(s,r));
            exit(-1);
          }*/
          //        id[l] = r;
        } else if(s->id[r] == 0) {
          if(l-r > -128) {
            s->id[r] = (signed char) l-r;
          } else {
            s->id[r] = (signed char)-128;
            s->idvtab = ALLOCMEMORY(space, s->idvtab, PairSint, s->idvcnt+1);
            s->idvtab[s->idvcnt].a = r;
            s->idvtab[s->idvcnt].b = l-r;
            s->idvcnt = s->idvcnt+1;
          //  qsort(s->idvtab, s->idvcnt, sizeof(PairSint), cmp_PairSint_qsort);
          }

          /*if (id(s,r) != l-r) {
            printf("[%d,%d]=[%d,%d]", l, r, r+id(s,r), r);
            exit(-1);
          }*/
         //        id[r] = -l;
        } else {
          printf("an error occured in id computation\n");
          exit(-1);
        }

        enqueue(space, &q, l);
        enqueue(space, &q, r);
      }
    }

    wrapList(space, list, destructinterval); 
 
  }
   
  qsort(s->idvtab, s->idvcnt, sizeof(PairSint), cmp_PairSint_qsort);
  wrapqueue(space, &q);
  return;
}

/*------------------------------- traverselcp --------------------------------
 *    
 * @brief performs a traversal of the suffix array
 * @author Steve Hoffmann 
 *   
 */

Uint *
traverselcp(void *space, Suffixarray *s){
    
  int i,
       lb,
       llcp,
       llb,
       min1,
       min2,
       m;

  Uint *A;
  
  unsigned char isLeaf; 
  Stack *stack,
        *mstack;
  
  A = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+2);
  memset(A, 255, sizeof(int)*s->numofsuffixes+2);

  stack = ALLOCMEMORY(space, NULL, Stack, 1);
  mstack = ALLOCMEMORY(space, NULL, Stack, 1);
  
  initStack(space, stack, 100000);
  initStack(space, mstack, 100000);

  /*push lcp and lbound*/
  stackpush(space, stack, 0);
  stackpush(space, stack, 0);

  for(i=1; i < s->numofsuffixes; i++) {
    lb = i-1;
    min1 = s->numofsuffixes+1;
    min2 = s->numofsuffixes+1;
    
    isLeaf = 1;    
    while (lcp(s,i) < stacktop(stack)) {
        isLeaf = 0;
        llcp = stackpop(stack);
        llb = stackpop(stack);
        
        /*child interval is given by llcp-[llb, i-1]*/
        /*cycle children here*/
        min1 = s->suftab[i-1];
        while (!stackisempty(mstack) && llb <= stacktop(mstack)) {
            stackpop(mstack);
            m = stackpop(mstack);
            
            if (m < min1) {
                min2 = min1;
                min1 = m;
            } else 
                if (m < min2){
                min2 = m;
            } 
        }
        
        lb = llb;

        stackpush(space, mstack, min1);
        stackpush(space, mstack, lb);
        if (id(s, lb) + lb == i-1) {
        //        if (id[lb] == i-1) {
            A[min2+1] = lb;
        } else {
            A[min2+1] = i-1;
        }
    }

    if(isLeaf) {  
      stackpush(space, mstack, s->suftab[i-1]);
      stackpush(space, mstack, lb);
    }

    if(lcp(s,i) > stacktop(stack)){
        stackpush(space, stack, lb);
        stackpush(space, stack, lcp(s,i));
    }
  }
  
  destructStack(space, stack);
  destructStack(space, mstack);
  return A;
}



/*--------------------------------- topdown ----------------------------------
 *    
 * @brief performs a top down traversal on the tree represented 
 * by the suffix array
 * @author Steve Hoffmann 
 *   
 */
 
void
constructsuflinks (void *space, Suffixarray *s, Uint *succ) {
  Uint i, 
       a, 
       b, 
       l, 
       r,
       u,
       v,
       w,
       m,
       n,
       d, 
       lidx,
       *B;
  List* list;
  Stack *stack;
   
  stack = ALLOCMEMORY(space, NULL, Stack, 1);  
  initStack(space, stack, 100000);

  B = ALLOCMEMORY(space, NULL, Uint, s->maxlcp+1);

  a = 0;
  b = s->numofsuffixes-1;
  stackpush(space, stack, b);
  stackpush(space, stack, a);


  while(!stackisempty(stack)) {
    a = stackpop(stack);
    b = stackpop(stack);


    list = getChildintervals(space, s, a, b);  
     
    d = getlcpval(s, a, b);
        
        if (id(s, a)+a == b) {
            B[d] = a;
        } else {
            B[d] = b;
        }

    for (i=0; i < list->length; i++) {

      l = ((PairUint*)list->nodes[i].data)->a;
      r = ((PairUint*)list->nodes[i].data)->b;         

      if(l != r) {

        stackpush(space, stack, r);
        stackpush(space, stack, l);
      
      } else {
      
        lidx = s->suftab[l];

        if(succ[lidx] < abs(s->numofsuffixes) && 
            (succ[lidx] > 0 || id(s, succ[lidx]) < s->numofsuffixes)) {
          
          if (id(s, succ[lidx]) > 0) {
            u = succ[lidx];
            v = id(s, succ[lidx]) + succ[lidx];
/*            v = id[succ[lidx]];*/
          } else {
/*            u = -1*id[succ[lidx]];*/
            u = succ[lidx] + id(s, succ[lidx]);
            v = succ[lidx];
          }
       
          d = getlcpval(s, u, v);

          if (id(s, B[d-1]) > 0) {
            m = B[d-1];
            n = id(s, B[d-1])+B[d-1];
            /*n = id[B[d-1]];*/
          } else {
            /*m = -1*id[B[d-1]];*/
            m = B[d-1] + id(s, B[d-1]);
            n = B[d-1];
          }
          
     /*     printf("%d-[%d,%d] -> %d-[%d,%d]\n", getlcpval(s, u, v), u, v, getlcpval(s, m, n), m, n);     
          for(w=0; w < 6; w++) {
            printf("%c", s->suffixptr[s->suftab[u]][w]);
          }
          printf(" -> ");
          for(w=0; w < 5; w++) {
            printf("%c", s->suffixptr[s->suftab[m]][w]);
          }
          printf("\n");
*/
  /*        for(w=0; w < getlcpval(s,u,v)-1; w++) {
            if(s->suffixptr[s->suftab[u]][w+1] != s->suffixptr[s->suftab[m]][w]){
                printf("mist\n");
                exit(-1);
            }
          }*/
        }
      }
    }

    wrapList(space, list, destructinterval); 
  }
  
  destructStack(space, stack);
  return;
}


inline void
addinterval(void *space, List *list, Uint a, Uint b) {
  PairUint *range;
  PairUint *check;

  range = ALLOCMEMORY(space, NULL, PairUint, 1);
  range->a=a;
  range->b=b;

  if(list->length > 0){
    check = (PairUint*) list->nodes[list->lastNode].data;
    if(range->a < check->a) {
      printf("check->a: %d, check->b: %d\n", check->a, range->a);
      while(1);
    }
  }

  if (list->length == 0) {
    insertAfter(space, list, LISTNILVALUE, range);
  } else {
    insertAfter(space, list, list->lastNode, range);
  }

  return;
}


Uint
cmp_leftbound(Uint a, void *toSort, void *key, void *info) {

  PairUint *intervala;
  PairUint *intervalb = NULL;

  List *list;
  Uint *l;

  l = (Uint*) key;

  list = (List*) toSort;

  intervala = (PairUint*) list->nodes[a].data;

  if(a < list->length-1)  {
    intervalb = (PairUint*) list->nodes[a+1].data;
  } else {
    return 0;
  }

  if (intervala->a <  *l && intervalb->a > *l)  return 0;
  if (intervala->a == *l && intervalb->a > *l)  return 0;
  if (intervala->a == *l && intervalb->a == *l) return 2;
  if (intervala->a <  *l && intervalb->a < *l)  return 2;

  if (intervala->a > *l) return 1;
  return 0;
}

List*
getChildintervals(void *space, 
    Suffixarray *s, 
    Uint i, 
    Uint j) { 

  List     *list; 
  Uint     i1,
           i2;
  unsigned char child;

  /*ALERT -1*/
  child = (i > 0 || j < s->numofsuffixes-1);
  list = initList(space, 10);


  if(child) {
    if (i < s->chldtab[j].val && s->chldtab[j].val <=j) {
      i1 = s->chldtab[j].val;
    } else {
      i1 = s->chldtab[i].val;
    }
    addinterval(space, list, i, i1-1);

  } else {
    i1 = i;
  }

  while(isnextlIndex(s,i1) && !isupIndex(s,i1) 
      && s->chldtab[i1].val != 0) {
    i2 = s->chldtab[i1].val;
    addinterval(space, list, i1, i2-1);
    i1 = i2;
  }

  if(child) {
    addinterval(space, list, i1,j);
  }

  return list;
}


inline Uint
getfirstlindex(Suffixarray *s, Uint i, Uint j){
  Uint val=0; 

  if((i==0 && j == s->numofsuffixes-1) || i==j) return 0;

  if (j < s->numofsuffixes && isupIndex(s,j) 
      && i < s->chldtab[j].val && j >= s->chldtab[j].val) {
    val = s->chldtab[j].val;
  } else if (isdownIndex(s,i)){
    val = s->chldtab[i].val;
  }

  return val;
}

inline Uint
getlcpval(Suffixarray *s, Uint i, Uint j){
  Uint val=0;

  if((i==0 && j == s->numofsuffixes-1) || i==j) return 0;

  if (j < s->numofsuffixes && isupIndex(s,j) 
      && i < s->chldtab[j].val && j >= s->chldtab[j].val) {
    val = lcp(s, s->chldtab[j].val);
  } else if (isdownIndex(s,i)){
    val = lcp(s, s->chldtab[i].val);
  }

  return val;
}


  PairUint
getCharInterval(void *space,
    Suffixarray *s,
    Uint i,
    Uint j,
    Uint pos,
    char ch) 
{
  List *list;
  Uint lcp=0;
  PairUint lr;

  lr.a = 1;
  lr.b = 0;

  if(i==j) return lr;

  list = getChildintervals(space, s, i, j);
  lcp = getlcpval(s, i, j);

  for(i=0; i < list->length; i++) {

    if(s->suffixptr[s->suftab[((PairUint*)list->nodes[i].data)->a]][lcp] == ch) {
      lr.a = ((PairUint*)list->nodes[i].data)->a;       
      lr.b = ((PairUint*)list->nodes[i].data)->b;
      break;
    }
  }
  wrapList(space, list, destructinterval);
  return lr;
}



  PairUint 
findslinkinterval( void *space, 
    Suffixarray *s, 
    List **lists,
    Uint listno, 
    Uint elem) 
{
  Uint  link;
  PairUint ret,
           test;
  List  *list;
  Uint obj;

  link = 
    s->inv_suftab[s->suftab[(((PairUint*) lists[listno]->nodes[elem].data)->a)]+1]; 

  /*printf("found link %d\n", link); */

  list = lists[listno-1];
  obj = binarySearch_m(list, list->length, &link, cmp_leftbound, NULL); 

  if (obj < list->length) {
    /*     printf("found elem: %d\n", elem);*/

    ret.a = ((PairUint*) list->nodes[obj].data)->a;
    ret.b = ((PairUint*) list->nodes[obj].data)->b;

    if (ret.a > link) {
      obj--;
      ret.a = ((PairUint*) list->nodes[obj].data)->a;
      ret.b = ((PairUint*) list->nodes[obj].data)->b;
    }

    while (++obj < list->length) {
      test.a = ((PairUint*) list->nodes[obj].data)->a;
      test.b = ((PairUint*) list->nodes[obj].data)->b;
      if(test.a <= link) {
        ret.a = test.a;
        ret.b = test.b;

      } else {
        break;
      }
    }

    /*printf("found %d, %d \n",ret.a, ret.b);*/
  } else {
    printf("nottin' found, maid!\n");
    ret.a = 0;
    ret.b = 0;
  }

  return ret;
}

  PairUint 
linkinterval( void *space, 
    Suffixarray *s, 
    List **lists,
    Uint i, 
    Uint ll) 
{
  Uint  link;
  PairUint ret,
           test;
  List  *list;
  Uint obj;

  link = s->inv_suftab[s->suftab[i]+1]; 

  /*printf("found link %d\n", link); */

  list = lists[ll-1];
  obj = binarySearch_m(list, list->length, &link, cmp_leftbound, NULL); 

  if (obj < list->length) {
    /*     printf("found elem: %d\n", elem);*/

    ret.a = ((PairUint*) list->nodes[obj].data)->a;
    ret.b = ((PairUint*) list->nodes[obj].data)->b;

    if (ret.a > link) {
      printf("inconsistency encountered\n");
    }
    while (++obj < list->length) {
      test.a = ((PairUint*) list->nodes[obj].data)->a;
      test.b = ((PairUint*) list->nodes[obj].data)->b;
      if(test.a <= link) {
        ret.a = test.a;
        ret.b = test.b;

      } else {
        break;
      }
    }

    /*printf("found %d, %d \n",ret.a, ret.b);*/
  } else {
    printf("nottin' found, maid!\n");
    ret.a = 0;
    ret.b = 0;
  }

  return ret;
}

void
constructsuflinks(void *space, Suffixarray *s) {

  Uint   i,
         total=0,
         j,
         a,
         b,
         k,
         nooflists,
         lcp,
         pos;
  Stack  istack;
  Stack  jstack;

  List   *children,
         **lists;
  PairUint **data,
           slinkinterval;

  nooflists = maxlcp(s) +1;
  /*printf(" allocating %d list ptr (%d Mbytes)\n", nooflists, 
    (int)(sizeof(List*)*nooflists)/1000000);*/

  lists = ALLOCMEMORY(space, NULL, List*, nooflists);
  memset(lists, 0, sizeof(List*)*nooflists);

  initStack(space, &istack, 1000);
  initStack(space, &jstack, 1000);

  stackpush(space, &istack, 0);
  stackpush(space, &jstack, s->numofsuffixes-1);

  while(!stackisempty(&istack)) {
    i = stackpop(&istack);
    j = stackpop(&jstack);
    lcp = getlcpval(s, i, j);

    /*printf("adding list %d\n", lcp);*/
    if (lists[lcp] == NULL) {
      lists[lcp] = initList(space, 5);
    }

    addinterval(space, lists[lcp], i, j);

    /*printf("lcp: %d-[%d,%d]\n", lcp, i, j);*/
    children = getChildintervals(space, s, i, j);
    data = (PairUint**) dataList(space, children);

    for(k=children->length; k > 0; k--) {

      a = data[k-1]->a;
      b = data[k-1]->b;

      FREEMEMORY(space, data[k-1]);

      if(a != b) { 
        total++;
        stackpush(space, &istack, a);
        stackpush(space, &jstack, b);
      }
    }

    FREEMEMORY(space, data);
    wrapList(space, children, NULL);
  }

  destructStack(space, &istack);
  destructStack(space, &jstack);

  s->suflink_l = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes);
  s->suflink_r = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes);
  memset(s->suflink_l, 0, sizeof(Uint)*s->numofsuffixes);
  memset(s->suflink_r, 0, sizeof(Uint)*s->numofsuffixes);

  for(i=1; i < nooflists-1; i++) {

    if(lists[i] !=  NULL && lists[i-1] !=NULL) {

      for(j=0; j < lists[i]->length; j++) {
        slinkinterval = findslinkinterval(space, s, lists, i, j);
        pos = getfirstlindex(s, ((PairUint*)lists[i]->nodes[j].data)->a, ((PairUint*)lists[i]->nodes[j].data)->b);
        s->suflink_l[pos]=slinkinterval.a;
        s->suflink_r[pos]=slinkinterval.b;
      }
    }
    wrapList(space, lists[i-1], destructinterval);
  }

  FREEMEMORY(space, lists);
  return;
}


void
suflinks(void *space, Suffixarray *s) {

  Uint   i,
         j,
         k,
         l,
         r,
         lk=0,
         ll,
         nooflists,
         pos;
  Queue  q;
  PairUint **data,
           slink;
  List *children;
  List **lists;
  unsigned char start = 1;

  nooflists = maxlcp(s) +1; 
  lists = ALLOCMEMORY(space, NULL, List*, nooflists);
  s->suflink_l = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes);
  s->suflink_r = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes);

  memset(s->suflink_l, 0, sizeof(Uint)*s->numofsuffixes);
  memset(s->suflink_r, 0, sizeof(Uint)*s->numofsuffixes);
  memset(lists, 0, sizeof(List*)*nooflists);
  emptyqueue(space, &q, 1000);

  enqueue(space, &q, 0);
  enqueue(space, &q, s->numofsuffixes-1);

  lists[0] = initList(space, 5);
  addinterval(space, lists[0], 0, s->numofsuffixes-1);

  while(!queueisempty(&q)) {
    i = dequeue(&q);
    j = dequeue(&q);

    children = getChildintervals(space, s, i, j);
    data = (PairUint**) dataList(space, children);
    DEBUG("processing interval %d-[%d,%d]\n", getlcpval(s,i,j), i, j);
    for(k=0; k < children->length; k++) {

      l = data[k]->a;
      r = data[k]->b;
      FREEMEMORY(space, data[k]);
      printf("adding child %d-[%d,%d]\n", getlcpval(s,l,r), l, r);
      if(l != r) {
        ll = getlcpval(s,l,r);
        enqueue(space, &q, l);
        enqueue(space, &q, r);


        if (ll==3) printf("adding to list 3 for [%d,%d] child [%d,%d]\n", i, j, l, r);
        if(lists[ll] == NULL) { 
          lists[ll] = initList(space, 500);
        }
        addinterval(space, lists[ll], l, r); 
      }
    }

    FREEMEMORY(space, data);
    wrapList(space, children, NULL);

    if (!start) {
      lk = getlcpval(s,i,j);
      if (lists[lk-1] != NULL && lists[lk-1]->length > 0) {
        if(i > 10000) {
          data = (PairUint**) dataList(space, lists[lk]);
          DEBUG("huhu %d %d", lk, i);

          for(k=0; k < lists[lk]->length; k++) {
            printf("%d-[%d,%d]\n", getlcpval(s, data[k]->a, data[k]->b), data[k]->a, data[k]->b);
          }
          printf("\n");
        }
        slink = linkinterval(space, s, lists, i, lk);
        pos = getfirstlindex(s, i, j);
        s->suflink_l[pos] = slink.a;
        s->suflink_r[pos] = slink.b;
      }
    }

    if(lk >= 2 && lists[lk-2] != NULL) {
      /*        wrapList(space, lists[lk-2], destructinterval);
                lists[lk-2] = NULL;*/
    }

    start = 0;
  }

  FREEMEMORY(space, lists);
  wrapqueue(space, &q);
  return;
}






void
dumplcptab(Suffixarray *s) {

  Uint i;

  for(i=0; i < s->numofsuffixes; i++) {
    printf("i:%d lcp:%d\n", 
        i, lcp(s,i));
  }

  printf("\n");

}


void
dumpchildtab(Suffixarray *s) {
  Uint i;

  for(i=0; i < s->numofsuffixes; i++) {
    printf("i:%d up:%d, down:%d, nextlIndex:%d := %d\n", 
        i, isnextlIndex(s,i), isdownIndex(s,i), isnextlIndex(s,i), s->chldtab[i].val);
  }

  printf("\n");
}



