
/*
 * sort_test.c
 * test the sort implementation
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include "memory.h"
 #include "sort.h"

 int main (int argc, char** argv) {
	Uint *sorted,i, *sorted2;
	char **ptr;
	int toSort[]={235657,2,3,0,67,234,1234,567,346,1,6,34,45,234,45,22,3,1,7,4,346,238,214,64,78,120,333};
  	const char x[]="asimple#treewillbegenerated#fromthis#string$";
    Uint len = (Uint) strlen(x);
	
  	printf("01234567890123456789012345678901234567890\n");
	printf("%s\n", x);
	printf("%c\n", x[len-1]);
	
	sorted=quickSortMultikey(NULL, &x, len-1, compareMkstr, len-1, NULL); 
	for(i=0; i < len; i++) printf("[%d]:%d =%c%c\n",i,sorted[i], x[sorted[i]],  x[sorted[i]+1]);
	printf("----------------------------------------\n");	
	ptr = ALLOCMEMORY(NULL, NULL, char *, len);
	for(i=0; i < len; i++) ptr[i] = (char*) &x[i];
	sorted2=quickSortMultikey(NULL, ptr, len-1, compareMkstrptr, len-1, NULL); 
	for(i=0; i < len; i++) 
	  printf("[%d]:%d =%c%c (%d: = %c%c)\n",i, sorted2[i], *ptr[sorted2[i]],  x[sorted2[i]+1], sorted[i], *ptr[sorted[i]],  x[sorted[i]+1]);
	
	printf("----------------------------------------\n");	
	 


	
	for(i=0; i < 25; i++) printf("[%d]:%d\n",i,toSort[i]);
	sorted=quickSort(NULL, &toSort, 25, cmp_int, NULL);
	for(i=0; i < 25; i++) printf("[%d]:%d\n",i,sorted[i]);
	for(i=0; i < 25; i++) printf("[%d]:%d\n",i,toSort[sorted[i]]);
	printf("\n");

	return (EXIT_SUCCESS);
 }


