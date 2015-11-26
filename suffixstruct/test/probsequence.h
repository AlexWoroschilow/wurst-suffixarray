#ifndef INTSEQUENCE_H
#define INTSEQUENCE_H

/*
 * probsequence.h
 * declaration of a probability sequence
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 */

 #include "basic-types.h"
 
 typedef struct {
	Uint *sequence;
	Uint length;
	
 } IntSequence;

#endif
