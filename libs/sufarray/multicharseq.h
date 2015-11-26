 #ifndef MULTI_SEQ_H
 #define MULTI_SEQ_H

/*
 *	multiseq.h
 *  declarations for a datastructure containing
 *  multiple integer sequences
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/11/06 15:09:15 CET  
 *
 */

 #include "basic-types.h" 
 #include "charsequence.h"

 typedef struct {
	void *ref; 	
	BOOL refersFile;

 } SeqReference;


 typedef struct {
 	Uint numofsequences;
	Uint totallength;
	Uint *markpos;		/*markpos[i] is the position of a*/
						/*separator character between S_i and S_i+1*/
	char *sequences; 	/*array of concatenated sequences*/
    SeqReference *ref;  /*ref[i] points to the original sequence*/
						/*that starts at position markpos[i]*/
    char *map;
    Uint mapsize;

 } MultiCharSeq;


 void dumpMultiCharSeq (MultiCharSeq *);
 MultiCharSeq* concatCharSequences(void *, CharSequence **, Uint, char, char);
 void destructMultiCharSeq(void*, MultiCharSeq *);
 Uint getMultiCharSeqIndex(MultiCharSeq *, char *);
 Uint getMultiCharSeqRelPos(MultiCharSeq *, char *);

 #endif
