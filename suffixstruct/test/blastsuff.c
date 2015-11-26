
/*
 *  blastsuff.c
 *  snippets
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 03/25/07 23:55:27 CEST
 *  
 */

in ranksufmatch 

			sum=0;
			printf("found %d positive matches of length %d\n", occ[i-1].m, 
				((imbissinfo*)info)->substrlen);
			explambda = ((imbissinfo*)info)->K * 
				  ((imbissinfo*)info)->substrlen * explambda; 
			
			printf("explambda: %19.16e\n", explambda);	
			for(v=0; v < occ[i-1].m; v++) {	
				explambda1 = explambda;
			   for(w=0; w < v; w++) {
				    printf("*\n");
					explambda1 = explambda1 * explambda1;
				}
				printf("|=%19.16e\n", explambda1);
			   sum += explambda1/ (v > 0 ? fak(v): 1);
			   printf("sum= %19.16e fak(v)=%d\n", sum, fak(v));
			}
			
			sum *= 1-exp(-explambda);
			
			printf("sum: %19.16e\n", sum);
			 
