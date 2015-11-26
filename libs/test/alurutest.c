
/*
 *  alurutest.c
 *  test aluruSort.c
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/11/2007 07:05:27 PM CEST
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "memory.h"
#include "biofiles.h"
#include "fileio.h"
#include "stringutils.h"
#include "aluruSort.h"
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>

int main(int argc, char** argv) {
   
  Uint i;
  int *space = NULL;
  char *pattern=  "MISSISSIPPI$";
  char *classes;  
  Uint *sdistance; 
  Uint *A;
  Uint **R;
  Uint no, max;
  Alurubucket *bckts;
  vector_t** qdlist;

  R = ALLOCMEMORY(space, NULL, Uint *, 1);
  classes = classify(space, pattern, strlen(pattern));
  sdistance = Qdist(space, classes, strlen(pattern), 'S');
  A = getAluruArray(space, pattern, strlen(pattern), '\0');
 
  printf("%s\n%s\n", pattern, classes);

  for(i=0; i < strlen(pattern); i++) {
    printf("%d",sdistance[i]);
  }

  printf("\n");
 
  max = uarraymax (sdistance, strlen(pattern)); 
  bckts = getAluruBuckets(space, pattern, strlen(pattern), &no, R);    
  qdlist = getQdistList(space, bckts, no, sdistance, strlen(pattern));

  
  showAluruBuckets(bckts, R[0], no);
  showQDlist(qdlist, sdistance[max]+1);
sortAluruSubstrings(space, qdlist, sdistance[max]+1, bckts, no, R[0], classes, strlen(pattern), 'S');

  return EXIT_SUCCESS;
}


