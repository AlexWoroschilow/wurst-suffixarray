
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

/* ------------------------------ cmpCharSequence ----------------------------
 *    
 * function to compare CharSequences for mulitkey sort (sort.c)
 * 
 */

  Uint
cmpCharSequence (Uint a, Uint b, Uint depth, void *data, void *info)
{
  char **s = (char**) data;	

  /*quick fix to meet end of multiintsequence criterion*/
  if(s[b][0] == (char) 127) {
    if (s[a][depth] == (char) 127) {
      return 0;
    }
    return 1;
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
  Suffixarray *arr;

  mseq = concatCharSequences(space, s, len, (char)254, (char)127);
  numofsuffixes = (mseq->totallength - mseq->numofsequences)+1; 

  printf("allocating space for %d suffixes (%d bytes)\n", numofsuffixes, 2*numofsuffixes*sizeof(Uint));
  inv_suftab = ALLOCMEMORY(space, NULL, Uint , numofsuffixes);
  arr = ALLOCMEMORY(space, NULL, Suffixarray, 1);


  printf("constructing suftab.\n");
  suffixptr = getSuffixPtr(space, mseq, (char)254, (char)127);
  sorted = quickSortMultikey (space, suffixptr, numofsuffixes, 
      cmpCharSequence, numofsuffixes-1, NULL);

  printf("constructing inv_suftab.\n");
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


/*------------------------------ destructSufArr ------------------------------
 *    
 * destruct a suffix array.
 * 
 */
  void
destructSufArr (void *space, Suffixarray *arr)
{
  FREEMEMORY(space, arr->suftab);
  if (arr->lcptab != NULL)
    FREEMEMORY(space, arr->lcptab);
  if (arr->inv_suftab != NULL)
    FREEMEMORY(space, arr->inv_suftab);
  if (arr->suffixptr != NULL)
    FREEMEMORY(space, arr->suffixptr);
  if (arr->seq != NULL)
    destructMultiCharSeq(space, arr->seq);
  FREEMEMORY(space, arr);

  return ;
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
  Uint s,t;
  int l=0;	

  arr->lcptab = ALLOCMEMORY(space, NULL, Uint, arr->numofsuffixes);


  for(i=0; i < arr->numofsuffixes; i++) {
    j = arr->inv_suftab[i];

    if (j > 0) {
      k = arr->suftab[j-1];
      s = arr->suffixptr[k]-arr->seq->sequences;
      t = arr->suffixptr[i]-arr->seq->sequences;

      l=l-1;
      if (l < 0) l=0;

      while (arr->seq->sequences[s+l] == arr->seq->sequences[t+l]){ 
        l++;
      }
      arr->lcptab[j] = l;
    }
  }
  arr->lcptab[0]=0;
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
        arr->lcptab[i],
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
    if (arr->lcptab[i] > 0) {
      s = arr->suffixptr[arr->suftab[i-1]]-arr->seq->sequences;
      t = arr->suffixptr[arr->suftab[i]]-arr->seq->sequences;
      printf("lcp of suffix %d and %d has length %d\t:\n", i-1, i, arr->lcptab[i]);
      for(j=0; j <= arr->lcptab[i]; j++) printf(" %d ", arr->seq->sequences[s+j]);
      printf("\n");
      for(j=0; j <= arr->lcptab[i]; j++) printf(" %d ", arr->seq->sequences[t+j]);
      printf("\n");
    }
  }
}

inline unsigned char
isnextlIndex(Suffixarray *s, Uint i) {
    return (s->lcptab[s->chldtab[i].val] == s->lcptab[i]);
}

inline unsigned char
isdownIndex(Suffixarray *s, Uint i) {
    return (s->lcptab[s->chldtab[i].val] > s->lcptab[i]);
}

inline unsigned char
isupIndex(Suffixarray *s, Uint i) {
    return (s->lcptab[i] > s->lcptab[i+1]);
}

void
constructchildtab(void *space, Suffixarray *s) {

  Uint i;
  int lastIndex = -1;
  Stack *stack;

  s->chldtab = ALLOCMEMORY(space, NULL, childtab, s->numofsuffixes);
  memset(s->chldtab, 0, s->numofsuffixes*sizeof(childtab));
  stack = ALLOCMEMORY(space, NULL, Stack, 1);
  initStack(space, stack, 100000);

  stackpush(space, stack, 1);

  for(i=1; i < s->numofsuffixes; i++) 
  {
    while(s->lcptab[i] < s->lcptab[stacktop(stack)]) {
      lastIndex = stackpop(stack);
      if(s->lcptab[i] <= s->lcptab[stacktop(stack)] && 
          s->lcptab[stacktop(stack)] != s->lcptab[lastIndex])
      {
        s->chldtab[stacktop(stack)].down = lastIndex;
        
        if (s->chldtab[stacktop(stack)].val != 0) printf("down conflict\n");
        s->chldtab[stacktop(stack)].val  = lastIndex;
      }
    }
    if (lastIndex != -1) {
      s->chldtab[i].up = lastIndex;
      
      if (s->chldtab[i-1].val != 0) printf("up conflict\n");
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
    while(s->lcptab[i] < s->lcptab[stacktop(stack)]) {
      stackpop(stack);
    }
    if (s->lcptab[i] == s->lcptab[stacktop(stack)]) {
      lastIndex = stackpop(stack);
      s->chldtab[lastIndex].nextlIndex = i;  
      s->chldtab[lastIndex].val = i;
    }
    stackpush(space, stack, i);
  }

  return;
}


void
addinterval(void *space, List *list, Uint a, Uint b) {
  PairUint *range;

  range = ALLOCMEMORY(space, NULL, PairUint, 1);
  range->a=a;
  range->b=b;
 
  if (list->length == 0) {
    insertAfter(space, list, LISTNILVALUE, range);
  } else {
    insertAfter(space, list, list->lastNode, range);
  }

  return;
}

void
destructinterval(void *space, void *data) {
  FREEMEMORY(space, (PairUint*) data);
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
    if (intervala->a == *l && intervalb->b > *l)  return 0;
    if (intervala->a == *l && intervalb->b == *l) return 2;
    if (intervala->a <  *l && intervalb->b < *l)  return 2;

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

  child = (i > 0 || j < s->numofsuffixes-1);
  list = initList(space, 10);


  if(child) {
    if (i < s->chldtab[j+1].up && s->chldtab[j+1].up <=j) {
      i1 = s->chldtab[j+1].up;
      if (!isupIndex(s,j) ||  s->chldtab[j].val != i1)
      printf("testing up: %d, testing val: %d=%d\n", isupIndex(s, j), s->chldtab[j].val, i1);
    } else {
      i1 = s->chldtab[i].down;

      if (!isdownIndex(s,i) ||  s->chldtab[i].val != i1)
      printf("testing down: %d, testing val: %d=%d\n", !isnextlIndex(s, i), s->chldtab[i].val, i1);
    }
    addinterval(space, list, i, i1-1);
  } else {
    i1 = i;
  }

  while(s->chldtab[i1].nextlIndex != 0) {
    i2 = s->chldtab[i1].nextlIndex;

    if (!isnextlIndex(s,i1) ||  s->chldtab[i1].val != i2)
    printf("testing nextlIndex: %d, testing val: %d=%d\n", isnextlIndex(s, i1), s->chldtab[i1].val, i2);
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

  if (i < s->chldtab[j+1].up && j >= s->chldtab[j+1].up)
    if(s->chldtab[j+1].up != s->chldtab[j].val || !isupIndex(s,j)) {
        printf("in getfirstlindex: is not an up index!\n");
    }
  
  if (i < s->chldtab[j+1].up && j >= s->chldtab[j+1].up)
    return s->chldtab[j+1].up;

  return s->chldtab[i].down;
}

inline Uint
getlcpval(Suffixarray *s, Uint i, Uint j){

  if (i == 0 && j == s->numofsuffixes-1) return 0;

   if (j+1 < s->numofsuffixes && 
      i < s->chldtab[j+1].up && j >= s->chldtab[j+1].up) {
        if(!isupIndex(s,j) || s->chldtab[j+1].up != s->chldtab[j].val){ 
          printf("not an upIndex!\n");
        }
   } else {
        
        if( (isdownIndex(s,i) && s->chldtab[i].down != s->chldtab[i].val)){
            if ((s->chldtab[i].down != s->chldtab[s->chldtab[i].val].val)) {
          printf("not an downIndex! isdown %d, isup %d, isnext %d, %d=%d\n", isdownIndex(s,i), isupIndex(s,i), isnextlIndex(s,i), s->lcptab[s->chldtab[i].val], s->lcptab[s->chldtab[i].down]); 
        }
        }
   }

  
  if (j+1 < s->numofsuffixes && 
      i < s->chldtab[j+1].up && j >= s->chldtab[j+1].up)
    return s->lcptab[s->chldtab[j+1].up];

  return s->lcptab[s->chldtab[i].down];
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


inline Uint
maxlcp(Suffixarray *s) {
  Uint i;
  Uint max = 0;

  for(i=0; i < s->numofsuffixes; i++) {
    if (s->lcptab[i] > max) max = s->lcptab[i];
  }
  return max;
}


PairUint 
findslinkinterval( void *space, 
                   Suffixarray *s, 
                   List **lists,
                   Uint listno, 
                   Uint elem) 
{
    Uint  link;
    PairUint ret;
    PairUint test;
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
          printf("encounterd suffix link inconsistency. continue ..\n");
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
      lists[lcp] = initList(space, 10);
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

  for(i=1; i < nooflists; i++) {
    if(lists[i] !=  NULL && lists[i-1] !=NULL) {
      for(j=0; j < lists[i]->length; j++) {
       /*printf("looking at interval [%d,%d], list %d\n", ((PairUint*)lists[i]->nodes[j].data)->a, ((PairUint*)lists[i]->nodes[j].data)->b, i);*/
        slinkinterval = findslinkinterval(space, s, lists, i, j);
        pos = getfirstlindex(s, ((PairUint*)lists[i]->nodes[j].data)->a, ((PairUint*)lists[i]->nodes[j].data)->b);
       /*printf("store at %d: [%d,%d]\n", pos, slinkinterval.a, slinkinterval.b);*/
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
dumplcptab(Suffixarray *s) {

  Uint i;

  for(i=0; i < s->numofsuffixes; i++) {
    printf("i:%d lcp:%d\n", 
        i, s->lcptab[i]);
  }

  printf("\n");

}


void
dumpchildtab(Suffixarray *s) {
  Uint i;

  for(i=0; i < s->numofsuffixes; i++) {
    printf("i:%d up:%d, down:%d, nextlIndex:%d\n", 
        i, s->chldtab[i].up, s->chldtab[i].down, s->chldtab[i].nextlIndex);
  }

  printf("\n");
}



