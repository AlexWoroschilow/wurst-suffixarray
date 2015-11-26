
/*
 *  wurstimbiss.c
 *  perform a boyer-moore exact string
 *  matching on probability sequences and score
 *  the matches w/ Wurst additionally
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/07/06 22:55:11 CET
 *  
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <math.h>
 #include <getopt.h>
 #include <sys/types.h>
 #include <sys/times.h>
 #include <time.h>
 #include <string.h>
 #include <ncurses.h>
 #include "wurstimbiss.h"
 #include "blaststat.h"
 #include "imbissblast.h"
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
 #include "depictseqs.h"
 #include "salami.h"
 #include "probseqscore.h"
 #include "usage.h"
 #include "dpalign.h"

/*WURST include*/ 
 #include "read_seq_i.h"
 #include "score_mat_i.h"
 #include "prob_vec_i.h"
 #include "prob_vec.h"
 #include "coord_i.h"
 #include "pair_set.h"
 #include "pair_set_i.h"
 #include "score_probvec.h"
 #include "coord.h"
 #include "align_i.h"
 #include "matrix.h"
 #include "model.h"
 #include "cmp_dmat_i.h"
 #include "altscores.h"
 #include "lsqf.h"
 #include "sort.h"

 static int verbose_flag;
 static struct option long_options[] =
 {
	/* These options set a flag. */
	{"verbose", no_argument,       &verbose_flag, 1},
	{"brief",   no_argument,       &verbose_flag, 0},	
	{"database",    required_argument, 0, 'f'},
	{"query",   required_argument, 0, 'i'},
	{"substrlen", required_argument, 0, 'l'},
	{"cutoff", required_argument, 0, 'c'},
	{"depictsw", no_argument, 0, 'd'},
	{"batchfile", required_argument, 0, 'b'},
	{"percent", required_argument, 0, 'p'},
	{"subfile", required_argument, 0, 's'},
	
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


void
latexscores (void *space, Matchtype *m, IntSequence **s, Uint len, Uint match, 
		   void *info)
{ 
  int sw;
  double explambda, E; 
  FILE* fp;
  stringset_t *query;

  	query = tokensToStringset(space, "/.", 
		s[m->id]->url, strlen(s[m->id]->url));

	fp=fopen("/Users/steve/coding/perl/decodepdb90/1s3zA.txt","a+");
	fprintf(fp, "%s\n", query->strings[query->noofstrings-1].str);
	fclose(fp);	
	
	explambda = exp(- ((imbissinfo*)info)->lambda * m->blast );
	E =  ((imbissinfo*)info)->K*3000000*12*explambda;	


  sw = floor(m->swscore);
  printf("%s & %s & %d & %.2f &",s[m->id]->url, s[m->id]->description, sw, 
	  	log(E));
  latexWurstAlignment(space, m, s, len, ((imbissinfo*)info)->query);
  printf("\\\\ \n");
  
  destructStringset(space, query);
}




void
allscores (void *space, Matchtype *m, IntSequence **s, Uint len, Uint match, 
		   void *info)
{
  char *pic;
  double explambda, E;
  FILE* fp;
  char *bin;
  float rmsd;
  stringset_t *query, *name;

  name = tokensToStringset(space, "/.", 
	  (((imbissinfo*)info)->query)->strings[0].str, 
	  (((imbissinfo*)info)->query)->strings[0].len );
 
  bin = ALLOCMEMORY(space, NULL, char, 100);
  
  strcpy(bin,"/Users/steve/coding/perl/decodepdb90/set/");
  strncat(bin, name->strings[name->noofstrings-2].str, 5);
  strncat(bin, ".lst", 4);

  query = tokensToStringset(space, "/.", 
		s[m->id]->url, strlen(s[m->id]->url));

  
  printf("[%d]: score: %f, count: %d\n", match, m->score, m->count); 
  printf("%d\t%s\t%d\t", m->id, s[m->id]->url, 
		 				  m->count);
  pic= depictSequence(space, len, 20, m->pos, m->count,'*');
  printf("[%s]\n", pic);
  printf("%s\n", s[m->id]->description);
  printf("gapless sw: %f\n", m->swscore);
  FREEMEMORY(space, pic);
 

  rmsd=doWurstAlignment(space, m, s, len, ((imbissinfo*)info)->query);

  fp = fopen(bin,"a+");
  fprintf(fp, "%s\t%f\t%d\n", query->strings[query->noofstrings-1].str,
	  rmsd, s[m->id]->length);
  fclose(fp);	

  FREEMEMORY(space, bin);
 
	printf("S:%f\n", m->blast);
	printf("lambda*S %19.16e\n", m->blast *((imbissinfo*)info)->lambda);
	explambda = exp(- ((imbissinfo*)info)->lambda * m->blast );
	printf("exp(-lambda*S): %19.16e\n", explambda);	
	E =  ((imbissinfo*)info)->K*3000000*12*explambda;	
	printf("E=Kmn * exp(-lambda*S): %19.16e\n", E);	
	printf("log(E): %f\n", log(E));	
	printf("1-exp(-E): %19.16e\n", 1-exp(-E)); 

	destructStringset(space, name);
	destructStringset(space, query);
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
	unsigned char depictsw=0;
	
	Uint i, j, noofvecs=0, nooffreqs=0, noofqueries=0;
	Uint substrlen = 10;
	Uint cutoff = 5;
	char *seq, *vec, *bin;
	
	void *space = NULL;	
	double *scores;
	char *path, *id;
	char *pveclistfile=NULL;
	char *alphabetfile=NULL;	
	char *inputfile=NULL;
	char *batchfile = NULL;
	char *subfile=NULL;



	imbissinfo *imbiss;

	stringset_t **fn, **freq, *queryurl, **queries=NULL;
	Suffixarray *arr;	
	IntSequence **sequences;
	IntSequence *input;
	FAlphabet *alphabet = NULL;		
    PairSint *matches=NULL;
	Uint percent=0;
	
	time_t startsuf, endsuf; 
	double difsuf, difmatch, difrank;	

#ifdef MEMMAN_H 	
	Spacetable spacetab;
	initmemoryblocks(&spacetab, 100000);
	space = &spacetab;
#endif
	
	while(1) 
	{
		c=getopt_long(argc, argv, "s:p:b:f:a:i:l:c:dv", long_options, &optindex);
		if (c==-1) break;
		
		switch(c) {
		  	case 'v':
			  verbose_flag=1;	
			  break;
		    case 'd':
			  	depictsw = 1;
				break;
			case 'f':
				pveclistfile = optarg;	
				break;
			case 'a':
				alphabetfile = optarg;	
				break;
			case 'i':
				inputfile = optarg;
				noofqueries = 1;
				break;
			case 'l':
				substrlen = atoi(optarg);
				break;
			case 'c':
			  	cutoff = atoi(optarg);
			    break;
			case 'b':
				batchfile = optarg;
				break;
			case 'p':
				percent = atoi(optarg);
				break;
			case 's':
				subfile = optarg;
				break;
			default:
				usage(argv[0]);
				exit (EXIT_FAILURE);
		}
	}
	if (pveclistfile==NULL || (inputfile == NULL && batchfile==NULL)) {
	usage(argv[0]);
		exit (EXIT_FAILURE);
	}
	
	/*read batch file*/
	if (batchfile) {
		queries = readcsv(space, batchfile, "", &noofqueries);	
	}

	/*read substitution matrix*/
	if (subfile) {
		freq=readcsv(space, subfile,",", &nooffreqs);
		scores = ALLOCMEMORY(space, NULL, double, ((nooffreqs-1)*(nooffreqs-1)) );
		for(i=1; i < nooffreqs; i++) {
			for(j=1; j < nooffreqs; j++) {
				if(strcmp(SETSTR(freq[i],j),"inf")==0){
					MATRIX2D(scores, nooffreqs-1, i, j)=0;
				}else{
					MATRIX2D(scores, nooffreqs-1, i, j)=atof(SETSTR(freq[i],j));
				}
			}
		}
	}
	
	/*read alphabet*/	
	if (alphabetfile != NULL) {
		alphabet = loadCSValphabet(space, alphabetfile);
		sortMapdomain(space, alphabet);
    }

	
	/*load sequence database*/
	fn=readcsv(space, pveclistfile, "", &noofvecs);
	sequences = ALLOCMEMORY(space, NULL, IntSequence *, noofvecs);
	for(i=0; i<noofvecs; i++) 
	{		  
		sequences[i] = loadSequence(space, SETSTR(fn[i],0));		
	}

	for (i=0; i < noofvecs; i++) {	
	  	destructStringset(space, fn[i]);
	}
	FREEMEMORY(space, fn);
	
	
	/*construct the suffix array*/
	time (&startsuf);
	arr = constructSufArr(space, sequences, noofvecs, NULL);
    constructLcp(space, arr); 	
   	time (&endsuf);
	difsuf = difftime (endsuf, startsuf);


	/*do search*/
    for (i=0; i < noofqueries; i++) {
	  
	    /*get query form batchfile*/
	  	if (queries) {
			inputfile = SETSTR(queries[i],0);
		}	
		
		input = loadSequence(space, inputfile);
		printf(">IMBISS order delivered\n");	
		seq = printSequence (space, input, 60); 
		printf("%s\n",seq);
		FREEMEMORY(space, seq);
 
		
		if (percent != 0) {
			substrlen = ((double)((double)input->length/100)*(double) percent);
			if (substrlen < 5) substrlen = 5;
		}
		
		time (&startsuf);
		matches=sufSubstring(space, arr, input->sequence, input->length, substrlen);	 
		time (&endsuf);
		difmatch = difftime (endsuf, startsuf);

    	queryurl = tokensToStringset(space, "/.", inputfile, strlen(inputfile));
		/*path = "/Users/steve/salami_lib/pvecs/";
		id = queryurl->strings[queryurl->noofstrings-2].str;
		printf("len: %d\n", strlen(id));	
		
		bin = attachpath(space, id, strlen(id), path, strlen(path), ".bin", 4);	
		vec = attachpath(space, id, strlen(id), path, strlen(path), ".vec", 4);
		
		printf("%s\n", bin);
		printf("%s\n", vec);
		*/
		bin = ALLOCMEMORY(space, NULL, char, 40);
		vec = ALLOCMEMORY(space, NULL, char, 40);
	
		strcpy(bin,"/Users/steve/salami_lib/pvecs/");
		strncat(bin, queryurl->strings[queryurl->noofstrings-2].str, 5);
		strcat(bin, ".bin");
		strcpy(vec, "/Users/steve/salami_lib/pvecs/");
		strncat(vec, queryurl->strings[queryurl->noofstrings-2].str, 5);
		strcat(vec, ".vec");
	
		
		destructStringset(space, queryurl);
		
		queryurl = initStringset(space);
		addString(space, queryurl, bin, 39);
		addString(space, queryurl, vec, 39);

		imbiss = ALLOCMEMORY(space, NULL, imbissinfo, 1);
		getimbissblast(space, input, sequences, noofvecs, alphabet, imbiss);
		imbiss->query = queryurl;
		imbiss->substrlen = substrlen;
		imbiss->alphabet = alphabet;
		if (subfile) {
			imbiss->sub = createsubmatrix(scores, imbiss->score, nooffreqs-1);
		}
		time (&startsuf);
		rankSufmatchFct(space, arr, matches, input->length-substrlen, cutoff, 1000, substrlen, 
		sequences, allscores, input, imbiss, scores, depictsw);	 
		time (&endsuf);
		difrank = difftime (endsuf, startsuf);
	
		printf ("Building  the suffixtree has taken %f seconds.\n", difsuf);
		printf ("Match the suffixtree has taken %f seconds.\n", difmatch);
    	printf ("Rank  the suffixtree has taken %f seconds.\n", difrank);
	
		
		destructStringset(space, queryurl);
		destructSequence(space, input);
		FREEMEMORY(space, imbiss->score);
		FREEMEMORY(space, matches);		
	}
	
	for (i=0; i < noofvecs; i++) {
		destructSequence(space, sequences[i]);
	}
    FREEMEMORY(space, sequences);	
	destructSufArr(space, arr);
	
#ifdef MEMMAN_H
	activeblocks(space);
#endif

	
	printf("Goodbye.\n");	
	return EXIT_SUCCESS;
}

