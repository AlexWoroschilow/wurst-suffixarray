
/*
 * testboyermoore.c
 * test implementation of exact string matching
 * Boyer-Moore style for arbitrary alphabets
 *
 * @author Steve Hoffmann
 * @date Sat 3 Nov 2006
 *
 */

 #include <stdlib.h>
 #include <string.h>
 #include "boyermoore.h"
 #include "basic-types.h"
 #include "mathematics.h"

 int main(int argc, char **argv) {
	vector_t *res;
               /*01234567890123456789012345678901234567890*/	
	char *tmpl ="guthardguthardguthardisgutundgotthardauch";
	char *patt = "gott";
	
	res = boyerMoore(NULL, NULL, 256, 
					 tmpl, strlen(tmpl), patt, strlen(patt), getASCIIcode);
	dumpVector(res);
	
	return EXIT_SUCCESS;
 }
 
