/*
 * mathematics_test.c
 * to test mathematics library
 *
 * @author Steve Hoffmann
 * @date Wed 22 Nov 2006 
 * 
 */

#include <stdlib.h>
#include <stdio.h>
#include "mathematics.h"
#include "cantor.h"

 void scrabbleMatrix_int(int *M, int m, int n) {
	int i,j;
	
	
	for (i=0; i < m; i++) {
		for (j=0; j < n; j++) {
			MATRIX2D(M,n,i,j)=(i*n)+j;	
		}
	}
 }

 void scrabbleMatrix3D_int(int *M, int m, int n, int l) {
	int i,j,k;
	
	
	for (i=0; i < m; i++) {
		for (j=0; j < n; j++) {
			for (k=0; k < l; k++) {
				MATRIX3D(M,m,l,i,j,k)=(i*n)+j;	
			}
		}
	}
 }


 int main(int argc, char **argv) {
   Uint i;
   int *matrix;
   vector_t vector, *v;
  
   printf("Hi there - testing mathematics.c now:\n");
  
   /*init 2D int matrix, scrabble and dump*/
   matrix = INITMATRIX2D(NULL, 100, 20, sizeof(int));   
   scrabbleMatrix_int(matrix, 100, 20);
   dumpMatrix_int(matrix, 100, 20);
   free(matrix);
   
   /*init 3D int matrix, scrabble and dump*/
   matrix= INITMATRIX3D(NULL,3,3,10,sizeof(int));
   scrabbleMatrix3D_int(matrix,3,3,10);
   MATRIX3D(matrix,3,3,0,0,2)=99; 
   dumpMatrix3D_int(matrix, 3, 3, 10);
   free(matrix);
 
   /*cantorize and de-cantorize a pair of numbers*/
   vector.length=0;
   vector.elements=NULL;
   APPENDVEC(NULL, &vector, 235);
   APPENDVEC(NULL, &vector, 3);
   APPENDVEC(NULL, &vector, 4);
   APPENDVEC(NULL, &vector, 10);
   
   i=codeCantor(&vector);
   printf("cantorized vector:\n");
   dumpVector(&vector);
   printf("to number: %d\n", i);
   
   v = decodeCantor(NULL, i, 4);
   printf("decoding routine returned:\n");
   dumpVector(v);
  
  
   /*init a vector, append some values and dump*/
   vector.elements = NULL;   
   vector.length = 0;
   for(i=0; i < 100; i++) { 
   		appendvector(NULL, &vector, i+1);
   }
 
   dumpVector(&vector);	
 
   /*test vector macros*/
   printf("vector ... \n");
   dumpVector(v);
   REVERSEVEC(0,LENGTHVEC(v)-1,v);
   printf("reversed...\n");
   dumpVector(v);
   SWAPVEC(0,1,v);
   SWAPVEC(1,2,v);
   dumpVector(v); 
 
   printf("doing the permutation now\n"); 
   /*generate several permutations and dump*/
   for(i=0; i < 23; i++) {
   	nextPermutation(v);
   	dumpVector(v);
   }
   
   return EXIT_SUCCESS;
 }

