
/*
 *  sufprobseq.c
 *  perform a boyer-moore exact string
 *  matching on probability sequences
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/07/06 22:55:11 CET
 *  
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <getopt.h>
 #include <sys/types.h>
 #include <sys/times.h>
 #include <time.h>
 #include "memory.h"
 #include "prob_vec.h"
 #include "prob_vec_i.h"
 #include "falphabet.h"
 #include "intsequence.h"
 #include "mathematics.h"
 #include "stringutils.h"
 #include "fileio.h"
 #include "createalphabet.h"
 #include "sufarray.h"
 #include "sufmatch.h"
 #include "mm.h"
 #include "encodeprobvec.h"
 
 const char *prog_bug="";
 const char *null_point="";

 static int verbose_flag;
 static struct option long_options[] =
 {
	/* These options set a flag. */
	{"verbose", no_argument,       &verbose_flag, 1},
	{"brief",   no_argument,       &verbose_flag, 0},	
	{"file",    required_argument, 0, 'f'},
	{0, 0, 0, 0}
 };


/*------------------------------- getProbChar --------------------------------
 *    
 * returns the position of character from a probability sequence over the 
 * given alphabet of size asize
 * 
 */
 
Uint
getProbChar (void *alphabet, Uint aszize, void *p, Uint pos)
{
  	/*at this point of the development alphabet is a stub*/
  	/*FAlphabet *a = (FAlphabet *)alphabet;*/
	IntSequence *seq = (IntSequence *)p;
	Uint ch = seq->sequence[pos];
	
	return ch;
}

/*---------------------------------- usage -----------------------------------
 *    
 * the usage function
 * 
 */
 
void
usage (char *name)
{
	printf("usage: %s -f <filename>\n", name);
  	return ;
}

/*----------------------------------- main -----------------------------------
 *    
 *  main function
 * 
 */
 
int
main (int argc, char** argv)
{	
 	Sint optindex, c;
	Uint i, noofvecs=0;
#ifdef MEMMAN_H	
	Spacetable spacetab;
#endif	  
	void *space = NULL;	
	char *pveclistfile;
	char *alphabetfile=NULL;	
	stringset_t **fn;
	Suffixarray *arr;	
	IntSequence **sequences;
	FAlphabet *alphabet = NULL;
	/*
	Uint j;
	Uint a[] = {1,1,1}  ;
	Uint b[] = {138,39,148,65,28,3,480,367,476,685,102,523,207,165,645,126,34,34,126};*/
    Uint immunoglobulin[] = {
396,  330,  932,  202,  202,  212,  212,  212,  212,   78,   54,   27,   52,
47,   5,  331,  328,  346,  346,  330,   62,   54,   27, 1022,   77,   92,
775,   66,   23,  493,  634,  395,   2,  328, 3,  551,  816,   54,  27,
52,   47,   5,  454,   50,   74,   56,   25,  675,  7,  90,  242,   80,
203,  174,  169,  28,  22,  367,  476, 71,  179,  222,  172,   63,  174,
169,  516,  332,  479,  616,  473,  154,  157,   10,   11,   14,    5,  328,
328,   62,   54,   27,  132,   77,   92,   71,   66,   23,  740,  865,   20,
330  ,347  ,351  ,555   ,32  ,574  ,454 ,1114  ,543  ,114  ,77  ,33  ,71,
66   ,23   ,24   ,51 ,1211  ,135 ,182  ,260  ,231  ,106  ,330    ,8
	};
	
	PairSint r;	
    PairSint *I;
	time_t startsuf, endsuf;
    clock_t clockstart, clockend;
	/*struct tms detailed;*/
	
	double difsuf;
	double dif;
	  
#ifdef MEMMAN_H 	
	initmemoryblocks(&spacetab, 100000);
	space = &spacetab;
#endif
	
    while(1) 
	{
		c=getopt_long(argc, argv, "f:a:", long_options, &optindex);
		if (c==-1) break;
		
		switch(c) {
			case 'f':
				pveclistfile = optarg;	
				break;
			case 'a':
				alphabetfile = optarg;	
				break;
			default:
				usage(argv[0]);
				exit (EXIT_FAILURE);
		}
	}
	if (pveclistfile==NULL) {
		usage(argv[0]);
		exit (EXIT_FAILURE);
	}
	

	fn=readcsv(space, pveclistfile, "", &noofvecs);
	if (alphabetfile != NULL) {
		alphabet = loadCSValphabet(space, alphabetfile);
		sortMapdomain(space, alphabet);
    }

	sequences = ALLOCMEMORY(space, NULL, IntSequence *, noofvecs);
	for(i=0; i<noofvecs; i++) 
	{		  
		sequences[i] = loadSequence(space, SETSTR(fn[i],0));		
	}

	time (&startsuf);
	arr = constructSufArr(space, sequences, noofvecs, NULL);
    constructLcp(space, arr); 	
   	time (&endsuf);
	difsuf = difftime (endsuf, startsuf);
	
	
	/*dumplcps(arr);*/
	/*
	clockstart = times(&detailed);
	printf("starting clock at %ld\n", (long) clockstart);


	//for(i=0; i < 1000000; i++) r=mmsearch(arr, a, 3, 0, 0, arr->numofsuffixes-1);
	//printf("found [%d, %d]\n", r.a, r.b);
	
	clockend = times(&detailed);
	printf("halted clock at %ld\n", (long)clockend);
	
	for (i=0; i <= (r.b-r.a)+2; i++) {
		printf("sequence %d:", i);
	    for (j=0; j < 3; j++) {
			printf(" %d ", arr->suffixptr[arr->suftab[r.a+i-1]][j]);
		}
		printf("\n");
	}*/


	r=mmsearch(arr, immunoglobulin, 100, 0, 0, arr->numofsuffixes-1);
    printf("immunoglobulin found at [%d, %d]\n", r.a, r.b);

	I=sufSubstring(space, arr, immunoglobulin, 115, 10);	
    reportSufmatch(arr, I, 115-10, 0, sequences);	
	rankSufmatch(space, arr, I, 115-10, 0, sequences);	
	/*destruct suffixarr*/

	/*destruct sequences*/
	  
    dif = ((double) clockend-clockstart)/(CLOCKS_PER_SEC);
	printf ("Building  the suffixtree has taken %.5lf seconds.\n", difsuf);
	printf ("Searching 1000000 patterns suffixtree has taken %.6lf s (%ld clk_tcks, %d clk_per_sec).\n", dif, (long) clockend-clockstart, (int) CLOCKS_PER_SEC);
	
	FREEMEMORY(space, arr);
	FREEMEMORY(space, fn);
	
#ifdef MEMMAN_H
	activeblocks(space);
#endif
	
	printf("Goodbye.\n");	
	return EXIT_SUCCESS;
}

