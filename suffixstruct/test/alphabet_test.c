/*
 * alphabet_test.c
 * testing alphabet routiness
 *
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include "memory.h"
 #include "falphabet.h"
 #include "mathematics.h"
 #include "cantor.h"
 #include "createalphabet.h"

 int main(int argc, char** argv) {
	FAlphabet *alphabet;
	vector_t* v;
	Uint ch,i;
	
    alphabet = loadCSValphabet(NULL, "alphabet_t100.csv");
	sortMapdomain(NULL, alphabet);
	dumpAlphabet(alphabet);
	ch=lookupChar(alphabet, 190344);
	
	printf("55331 found at position %d, lookup shows ch:%d and map:%d\n", ch, alphabet->characters[ch], alphabet->mapdomain[ch]);
	v = decodeCantor(NULL, alphabet->mapdomain[ch], 2);
	for (i=0; i < LENGTHVEC(v); i++) {
		printf("[%d]:%d\n", i, VECTOR(v,i));
	}
	
	return EXIT_SUCCESS;
 }
 

