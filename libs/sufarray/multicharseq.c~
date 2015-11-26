
/*
 *  multiseq.c
 *  some functions to handle multiseqs (type char)
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/15/06 11:42:53 CET
 *  
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "charsequence.h"
 #include "vtprogressbar.h"
 #include "multicharseq.h"
 #include "sort.h"
 
/*---------------------------- concatCharSequences ----------------------------
 *    
 * concatenates CharSequences using a given Uint delimiter
 * and stores them in a MultiCharSeq container.
 * 
 */
 
MultiCharSeq *
concatCharSequences (void *space, CharSequence **s, Uint len, 
					char delim, char sentinel)
{
    char *buf=NULL;
    char *map = NULL;
    Uint i, j, k=0, 
		 totallength=0, 
		 *markpos;
	MultiCharSeq *mseq;

	mseq = ALLOCMEMORY(space, NULL, MultiCharSeq, 1);
	markpos = ALLOCMEMORY(space, NULL, Uint, len);
    map = ALLOCMEMORY(space, NULL, Uint, 257);
    memset(map, 0, 256);

	for(i=0; i < len; i++) {
	  
        totallength += (s[i]->length+1);
		buf = ALLOCMEMORY(space, buf, char, totallength+1);
		if (buf==NULL) fprintf(stderr, "allocation failed: exiting\n");

//		progressBarVT("sequences processed", len-1, i, 25);	
		for(j=0; j < s[i]->length; j++) {
			buf[k] = s[i]->sequence[j];
            map[(Uint)buf[k]]=buf[k];
            k++;
		}
		/*separate sequences or finalize*/
		if (i == (len-1)) {
		  buf[k] = sentinel;
          map[(Uint)buf[k]]=buf[k];
		  markpos[i] = k;
		  k++;
          buf[k]='\0';

		} else {
		  buf[k] = delim;
		  map[(Uint)buf[k]]=buf[k];
          markpos[i] = k;
		  k++;
		}
        
        /*FREEMEMORY(space, s[i]->sequence);*/
	}
	mseq->totallength = totallength;
	mseq->numofsequences = len;
	mseq->sequences = buf;
	mseq->markpos = markpos;

    for(i=0; i < 256; i++) {
        if(map[i]==0) {
           j=i+1;
           while(j<256 && map[j]==0) j++;
           if (j < 256) {
             map[i]=map[j];
             map[j]=0;
           } else {
             break;
           }
        }
    }

    map = ALLOCMEMORY(space, map, char, i+1);
    mseq->map = map;
    mseq->mapsize = i;


	return mseq;
}


/*----------------------------- destructMultiSeq -----------------------------
 *    
 * destructs a MultiSeq structure
 * 
 */

void
destructMultiCharSeq (void *space, MultiCharSeq *mseq)
{
    
	FREEMEMORY(space, mseq->sequences);
	FREEMEMORY(space, mseq->markpos);
	FREEMEMORY(space, mseq);
	return ;
}


/*------------------------------- cmp_markpos --------------------------------
 *    
 * compare function for getMultiSeqIndex
 * 
 */
 
Uint
cmp_markpos (Uint a, void *data, void *key, void *info)
{
    Uint *d = (Uint*) data;
	Uint *k = (Uint*) key;
	
	if (d[a] > *k) {
		if (a > 0) 
		{
			if (d[a-1] < *k) 
			{
				return 0;
			} 
			else 
			{
				return 1;
			}
		} 
		else 
		{
			return 0;
		}
	}
	
    if (d[a] < *k) return 2;

	return 0;
}

/*-------------------------- getMultiSeqDescription --------------------------
 *    
 * returns the index of a sequence in multiseq addressed by a pointer
 * 
 */
 
Uint
getMultiCharSeqIndex (MultiCharSeq *mseq, char *ptr)
{	
	Uint pos, i;
	
	pos = (ptr - mseq->sequences);
 
	i=binarySearch(mseq->markpos, mseq->numofsequences, &pos, 
			       cmp_markpos, NULL);

	return i;
}

/*---------------------------- getMultiSeqRelPos -----------------------------
 *    
 * returns the relative position of a pointer to mulitseq
 * with respect to the addressed sequence.
 * 
 */
 
Uint
getMultiCharSeqRelPos (MultiCharSeq *mseq, char *ptr)
{
	return 0;
}


/*------------------------------- dumpMultiSeq -------------------------------
 *    
 * dumps a multiseq to the screen
 * 
 */

void
dumpMultiCharSeq (MultiCharSeq *mseq)
{
  	Uint i;

	for(i=0; i < mseq->totallength; i++) {
		printf("%c-", mseq->sequences[i]);	
	}

	printf("\n");
	return ;
}

