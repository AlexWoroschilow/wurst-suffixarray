/*
 * fileio_test.c
 * testing file_io routines
 *
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include "memory.h"
 #include "fileio.h"
 #include "stringutils.h"

 int main(int argc, char** argv) {
	char* content;
	Uint contentlen,i,j;
	stringset_t *set, *set2, **csv;
    int space;	


	content = readfile(&space, "casessmall.rtxt", &contentlen);
	printf("read file of length: %d\n", contentlen);
	set = tokensToStringset(&space, "\n", content, contentlen);
	FREEMEMORY(&space, content);
	printf("file contained %d lines\n", set->noofstrings);
	
	for(i=0; i < set->noofstrings; i++) {
		
		set2 = tokensToStringset(&space, ",", set->strings[i].str, set->strings[i].len);
		for(j=0; j < set2->noofstrings; j++) {
			printf("%s,", set2->strings[j].str);
		}
		printf("\n");
		destructStringset(&space, set2);
	}
	destructStringset(&space, set);
		
	printf("reading csv ...\n");
	csv = readcsv(&space, "cases.rtxt", ",", &contentlen);
	
	for (i=0; i < contentlen; i++) {
		for(j=0; j < csv[i]->noofstrings; j++) {
			printf("%s",csv[i]->strings[j].str);
		}
		destructStringset(&space, csv[i]);
		printf("\n");
	}
	
	FREEMEMORY(&space, csv);
	
	return EXIT_SUCCESS;
 }
 

