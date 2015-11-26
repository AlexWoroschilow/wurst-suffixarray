
/*
 *  aligntest.c
 *  test for dpalign methods
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 03/19/07 02:15:23 CET
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "dpalign.h"

void testedist(void *space) {

	int i, j, *sub;
  	char *sb = "dbadad";
	char *sa = "ada";
	
	sub = INITMATRIX2D(space, 1500, 1500, sizeof(Uint));
	sub = memset(sub, 0, sizeof(Uint)*(1500*1500));
	for(i=0; i < 300; i++) {
		for(j=0; j < 300; j++) {
			MATRIX2D(sub, 1500, i, j) = 1;
		}
	}
	
	printf("the distance is %d\n", edist(space, sa, 3, sb, 6, 1, sub, 1500));
    FREEMEMORY(space, sub);
}


void swtest(void *space)
	int swscores[2] = {3,-2};
	int *swres;
  	char *sb = "dbadad";
	char *sa = "ada";
	
	swres = swmatrix(space, sb, 6, sa, 3, -10, constscr, swscores); 
	dumpMatrix_int(swres, 6+1, 3+1);
	printf("max: %d\n", swres[arraymax(swres, (6+1)*(3+1))]);
	FREEMEMORY(space, swres);
}


