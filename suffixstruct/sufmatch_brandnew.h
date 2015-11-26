#ifndef SUF_MATCH
#define SUF_MATCH

/*
 *
 *	sufmatch.h
 *  declarations for matching functions in suffixarrays
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/19/06 15:16:31 CET  
 *
 */

 #include "basic-types.h"
 #include "intsequence.h"
 #include "sufmatch.h"
 #include "sufarray.h" 

typedef struct {
	Uint id;
	Uint count;
	Uint *pos;
	Uint *org;
	Uint m;
	float score;
	float swscore;
	double blast;
} Matchtype;

PairSint* sufSubstring (void *, Suffixarray *, Uint *, Uint, Uint);
void reportSufmatch    (Suffixarray *, PairSint *, Uint, Uint, 
						IntSequence **);


void
rankSufmatch ( 	void *space, 
				Suffixarray *a, 
				PairSint *matches, 
				Uint len, 
				Uint threshold, 
				Uint K, 
				Uint S, 
				IntSequence **s,
				double (*fltr)(		void *, 
				  					Matchtype *, 
									IntSequence *, 
									IntSequence *, 
	  				        		Uint *, 
									Uint, 
									Uint, 
									void *),
				Matchtype* (*sel)(	void *, 
				  					Matchtype *, 
									Uint,
									IntSequence *,
									IntSequence **,
									void *),
				void (*handler)(	void *, 
				  					Matchtype *, 
									IntSequence **, 
									Uint, 
									Uint, 
									void *), 	
				IntSequence *queryseq, 
				void *info, 
				double *scores, 
				unsigned char depictsw
);



Matchtype*
selectBlastSWconst(void *space, Matchtype *m, Uint k, 
	IntSequence *a, IntSequence **s, void *info);
 
Matchtype*
selectSW (void *space, Matchtype *m, Uint k, IntSequence *a, 
		IntSequence **s, void* info);


Matchtype*
selectBlast (void *space, Matchtype *m, Uint k, IntSequence *a, 
		IntSequence **s, void* info);

double
swconstfilter(void *space, Matchtype *m, IntSequence *a, IntSequence *b,
					Uint *ptr, Uint len, Uint pos, void *info); 

double
scorefilter(void *space, Matchtype *m, IntSequence *a, IntSequence *b,
					Uint *ptr, Uint len, Uint pos, void *info);


/*void rankSufmatch      (void *, Suffixarray *, PairSint *, Uint, Uint, 
						IntSequence **);
*/
						
void rankSufmatchList  (void *, Suffixarray *, PairSint *, Uint, Uint, 
						IntSequence **);
void rankSufmatchFct   (void *, Suffixarray *, PairSint *, Uint, Uint,
						Uint, Uint, IntSequence **, 
						void (*handler)(void *, Matchtype *, IntSequence **, 
						  				Uint, Uint, void*),
						IntSequence *, void *, double*, unsigned char);
#endif

