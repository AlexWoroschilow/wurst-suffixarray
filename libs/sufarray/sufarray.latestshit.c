
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

inline PairUint
getSuflink(Suffixarray *s, Uint i, Uint j) {
  Uint slidx;
  PairUint link;
  int a, b;

  slidx = getfirstlindex(s, i, j);

  if (id(s, s->suflink[slidx]) > 0) {
    a = s->suflink[slidx];
    b = id(s, s->suflink[slidx])+s->suflink[slidx];
  } else {
    a = s->suflink[slidx]+id(s, s->suflink[slidx]);
    b = s->suflink[slidx];
  }

  link.a = a;
  link.b = b;
  return link;
}

/*---------------------------- constructchildtab -----------------------------
 *    
 * @brief performs bottom-up traversals to construct childtable 
 * @author Steve Hoffmann 
 *   
 */
 
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



/*------------------------------- computeId ----------------------------------
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
          }

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
          } 
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


/*------------------------------- getsuffsucc --------------------------------
 *    
 * @brief performs a bottom-up traversal of the suffix array collecting
 *        cause(p)-successors for suffix link construction
 * @author Steve Hoffmann 
 *   
 */

Uint *
getsufsucc(void *space, Suffixarray *s){
    
  int  i,
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


  stackpush(space, mstack, s->suftab[0]);
  stackpush(space, mstack, 0);
  
  for(i=1; i < s->numofsuffixes; i++) {
    lb = i-1;
    
    stackpush(space, mstack, s->suftab[i-1]);
    stackpush(space, mstack, i-1);
    
    isLeaf = 1;    
    while (lcp(s,i) < stacktop(stack)) {
        isLeaf = 0;
        llcp = stackpop(stack);
        llb = stackpop(stack);
        
        /*child interval is given by llcp-[llb, i-1]*/
        /*cycle children here*/

     /* min1 = (s->suftab[i-1] < s->suftab[llb]) ? s->suftab[i-1] : s->suftab[llb];
        min2 = (s->suftab[i-1] < s->suftab[llb]) ? s->suftab[llb] : s->suftab[i-1];
     */   
        
        
        min1 = s->numofsuffixes+1;
        min2 = s->numofsuffixes+1;
        while (!stackisempty(mstack) && llb <= stacktop(mstack) && 
                stacktop(mstack) <= i-1) {
            stackpop(mstack);
            m = stackpop(mstack);
            
            if (m < min1) {
                min2 = min1;
                min1 = m;
            } else {
              if (m < min2 && m != min1)
                min2 = m;
            } 
        }        
        lb = llb;

        printf("interval %d-[%d,%d] - min1: %d, min2 %d; {%d,%d}\n", llcp, llb, i-1, min1, min2, s->inv_suftab[s->suftab[llb]+1],s->inv_suftab[s->suftab[i-1]+1]); 

        stackpush(space, mstack, min1);
        stackpush(space, mstack, lb);
        if (id(s, lb) + lb == i-1) {
            A[min2+1] = lb;
        } else {
            A[min2+1] = i-1;
        }
    }

/*    if(isLeaf) {  
      
      printf("singleton [%d,%d] - suf: %d\n", lb, i-1, s->suftab[i-1]);
      stackpush(space, mstack, s->suftab[i-1]);
      stackpush(space, mstack, lb);
    } else {
      stackpush(space, mstack, s->suftab[i-1]);
      stackpush(space, mstack, lb);

    } */

    if(lcp(s,i) > stacktop(stack)){
        stackpush(space, stack, lb);
        stackpush(space, stack, lcp(s,i));
    }
  }
  
  destructStack(space, stack);
  destructStack(space, mstack);
  return A;
}


/*---------------------------- constructsuflinks ----------------------------
 *    
 * @brief performs a top down traversal on the tree represented 
 * by the suffix array
 * @author Steve Hoffmann 
 *   
 */

void
constructsuflinks (void *space, Suffixarray *s, Uint *succ) {
  Uint i, 
       a, b, 
       l, r;
  int  u, v,
       m, n;
  Uint d, 
       slidx,
       lidx,
       *B; 
  List *list;
  Stack *stack; 
  PairUint suflink;
  
  stack = ALLOCMEMORY(space, NULL, Stack, 1);  
  initStack(space, stack, 100000);
  s->suflink = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+1);

  s->suflink_l = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+1);
  s->suflink_r = ALLOCMEMORY(space, NULL, Uint, s->numofsuffixes+1);
  
  memset(s->suflink_l, 0, s->numofsuffixes-1*sizeof(Uint));
  memset(s->suflink_r, 0, s->numofsuffixes-1*sizeof(Uint));
  
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
    } else if (a == b+id(s, b)){
      B[d] = b;
    } else {
        printf("wrong id\n");
    }

    printf("[%d,%d] has the following children: \n", a, b);
    for (i=0; i < list->length; i++) {

      l = ((PairUint*)list->nodes[i].data)->a;
      r = ((PairUint*)list->nodes[i].data)->b;         
      printf("-- [%d,%d]\n", l, r);

      if(l != r) {
         
        stackpush(space, stack, r);
        stackpush(space, stack, l);

      } else {

        lidx = s->suftab[l];
     
        printf("checking leaf: [%d,%d] lidx:%d, succ[lidx]:%d\n", l, l, s->suftab[l], succ[s->suftab[l]]);
        if((succ[lidx] > 0 && succ[lidx] < s->numofsuffixes-1) &&
            (abs(id(s, succ[lidx])) < s->numofsuffixes-1)) {  

          if (id(s, succ[lidx]) > 0) {
            u = succ[lidx];
            v = id(s, succ[lidx]) + succ[lidx];
          } else {
            u = succ[lidx] + id(s, succ[lidx]);
            v = succ[lidx];
          }

          d = getlcpval(s, u, v);

          if (id(s, B[d-1]) > 0) {
            m = B[d-1];
            n = id(s, B[d-1])+B[d-1];
     //       printf("pos [%d,%d]\n", m, n);
          } else {
            m = B[d-1] + id(s, B[d-1]);
            n = B[d-1];
     //       printf("neg [%d,%d]\n", m, n);
          }

          slidx = getfirstlindex(s, u, v);
          s->suflink[slidx] = B[d-1];
          
          s->suflink_l[slidx] = m;
          s->suflink_r[slidx] = n;
          suflink = getSuflink(s, u, v);

           
         printf("linking [%d,%d] -> [%d,%d] {%d,%d}\n", 
             u,v,m,n, s->inv_suftab[s->suftab[u]+1],s->inv_suftab[s->suftab[v]+1]); 

         #ifdef SUFLINKDEBUG
          { int w;
            for(w=0; w < getlcpval(s,u,v)-1; w++) {
              if(s->suffixptr[s->suftab[u]][w+1] != s->suffixptr[s->suftab[m]][w] || getlcpval(s, u, v) != getlcpval(s, m, n)+1){
                printf("mist\n");
                exit(-1);
              }
               if(s->suffixptr[s->suftab[v]][w+1] != s->suffixptr[s->suftab[n]][w] || getlcpval(s, u, v) != getlcpval(s, m, n)+1){
                printf("mist\n");
                exit(-1);
              }

            }
          }
          #endif
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



