
/*
 *  biofiletest.c
 *  test for the biofile source
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/10/2007 02:50:57 PM CEST
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "biofiles.h"
#include "fileio.h"
#include "stringutils.h"
#include "charsequence.h"
#include "multicharseq.h"
#include "sufarray.h"
#include "gnuplot_i.h"
#include "mmchar.h"
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>

int main(int argc, char** argv) {
  char* content;
  Uint contentlen, i, j, k, l, id, lines=0;
  stringset_t *set, *set2, **csv, *que;
  CharSequence **s;
  Suffixarray *sarray;
  MultiCharSeq *mseq;
  PairSint d, *matches  = NULL;
  Uint totallength = 0;
  Uint wsize=10;
  Uint counter=0;
  Uint all=0;
  int *space = NULL;
  char *pattern=  "GGAAGAAAGCGTGGGGTTTG";
  char *pattern2= "TGATTAGTGATTAGTGATTA";
  char *pattern3= "ACAAACATAT";
  char *start;
  time_t startsuf, endsuf; 
  double difsuf;
  Uint noofchildren;
  List *list;
  PairUint **childinterval;
  gnuplot_ctrl *h;
  double *genome;
  
  //set = readfasta(&space, "HP26695.fasta");
  //csv = readcsv(&space, "HP12_GCTC.inserts", "", &lines); 
  /*s = ALLOCMEMORY(&space, NULL, CharSequence *, set->noofstrings);   


   for(i=0; i < set->noofstrings/2; i++) {
    totallength += set->strings[i].len; 

    s[i] = ALLOCMEMORY(&space, NULL, CharSequence, 1);
    s[i]->sequence = set->strings[i].str;
    s[i]->length = set->strings[i].len;
    /  printf("%s,", set->strings[i].str);
        printf("\n"); / 
  }*/



  s = ALLOCMEMORY(&space, NULL, CharSequence *, 1);
  s[0] = ALLOCMEMORY(&space, NULL, CharSequence, 1);
  s[0]->sequence = pattern3;
  s[0]->length = strlen(pattern3);


  genome = ALLOCMEMORY(&space, NULL, double, totallength);
  memset(genome, 0, sizeof(double)*totallength);


  time (&startsuf);
  sarray = constructSufArr(&space, s, 1/*set->noofstrings/2*/, NULL); 
  constructLcp(space, sarray);
  dumplcptab(sarray);
  constructchildtab(space, sarray);
  time (&endsuf);
  difsuf = difftime (endsuf, startsuf);

  printf("noofsuffixes: %d\n", sarray->numofsuffixes);

  dumpchildtab(sarray);
  dumpSufArr(sarray);
  
  list = getChildintervals(space, sarray, 0, 5);
  childinterval = (PairUint**) dataList(space, list);
  for(i=0; i < list->length; i++) {
    printf("[%d,%d]\n", childinterval[i]->a, childinterval[i]->b);
  }
  
  constructsuflinks(space, sarray);

  for(k=1; k < lines; k+=2) {
//    printf("searching %s\n", csv[k]->strings[0].str);

  
    if(csv[k]->strings[0].len > 8) {  
    if(wsize > csv[k]->strings[0].len) {
      d=mmsearch(sarray, csv[k]->strings[0].str, csv[k]->strings[0].len, 0, 0, sarray->numofsuffixes-1);
  //    printf("suffixes were found at positions (%d, %d)\n",d.a, d.b);
        for  (j=d.a; j <= d.b; j++) {
            genome[sarray->suftab[j]]++;
         }

    } else {
      for(l=0; l < csv[k]->strings[0].len-wsize; l++) { 
         d=mmsearch(sarray, &csv[k]->strings[0].str[l], wsize, 0, 0, sarray->numofsuffixes-1);
    //    printf("suffixes were found at positions (%d, %d)\n",d.a, d.b);
        
        
     for  (j=d.a; j <= d.b; j++) {
            genome[sarray->suftab[j]]++;

       /*   start = sarray->suffixptr[sarray->suftab[j]];
          printf("pattern was: %s\n", &csv[k]->strings[0].str[l]);
          printf("suffix found: ");
          for (i=0; i < wsize; i++) {
          printf("%c", start[i]);
          }
          printf("\n");
          id = getMultiCharSeqIndex(sarray->seq, sarray->suffixptr[sarray->suftab[j]]);	
          printf("found in sequence: %d\n", id); */
        } 
      }
    }
    if (d.a < d.b) counter++;
    all++;}
  }
  
  destructStringset(&space, set);
  writeY("out.xy", genome, totallength);

  /*h = gnuplot_init();
  gnuplot_setstyle(h, "points");
  
  gnuplot_cmd(h, "set title 'IMBISS - seed statistics' -28,0 font'Helvetica,15'");	
  gnuplot_cmd(h, "set label 'seed length: %d' at graph 0.05,0.95 font 'Helvetica, 12'", totallength);
  gnuplot_set_xlabel(h, "matches");
  gnuplot_set_ylabel(h, "position");	
  gnuplot_plot_x(h, genome, totallength, "position");
  */

  printf ("sliding windows of %d sequences (of %d) found\n", counter, all);
  printf ("Building  the suffixarray has taken %f seconds.\n", difsuf);
  printf ("Total length of suffixarray was %d\n", totallength);
  while(1);
  return EXIT_SUCCESS;
}

