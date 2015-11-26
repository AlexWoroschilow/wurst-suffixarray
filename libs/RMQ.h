
/*
 *
 *	RMQ.h
 *  declarations for range minimum queries
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 10/01/2007 06:54:12 PM CEST  
 *
 */


int
RMQuery (int node, int b, int e, int* M, int* range, int i, int j);

int*
RMQinit (int node, int b, int e, int* M, int* range, int n);

