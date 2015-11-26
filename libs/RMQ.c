
/*
 *  RMQ.c
 *  processing range minimum queries
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 10/01/2007 03:10:23 PM CEST
 *  
 */

#include "mathematics.h"



/*--------------------------------- RMQinit ----------------------------------
 *    
 * @brief initializes a segment tree for RM queries.
 * @author Steve Hoffmann 
 *   
 */

void
RMQinit (int node, int b, int e, int* M, int* range, int n) {
  if (b==e) {
    M[node] = b;
  } else {
    RMQinit(2*node, b, (b+e)/2, M, range, n);
    RMQinit(2*node+1, (b+e)/2 +1, e, M, range, n);
  }

  if(range[M[2*node]] <= range[M[2*node+1]]) {
    M[node] = M[2*node];
  } else {
    M[node] = M[2*node+1];
  }
}



/*--------------------------------- RMQuery ----------------------------------
 *    
 * @brief perform a range minimum query
 * @author Steve Hoffmann 
 *   
 */

int
RMQuery (int node, int b, int e, int* M, int* range, int i, int j) {
   
   int p1, p2;
   
   if(i > e || j < b) 
      return -1;

    if(b >= i && e <= j)
      return M[node];

    p1 = RMQuery(2*node, b, (b+e)/2, M, range, i, j);
    p2 = RMQuery(2*node +1, (b+e)/2 +1, e, M, range, i, j);

    if (p1 == -1)
      return M[node] = p2;
    if (p2 == -1)
      return M[node] = p1;
    if (range[p1] <= range[p2])
      return M[node] = p1;

    return M[node]=p2;
  }


