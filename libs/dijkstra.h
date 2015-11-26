
/*
 *
 *	dijkstra.h
 *  
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 09/15/2007 03:30:39 PM CEST  
 *
 */

#define INFINITY    (MAX_INT - 1)

typedef struct {
  int weight;
  int dest;
} DijkEdge;

typedef struct {
  DijkEdge* connections; /* Un array de arcos */
  int numconnect;
  int distance;
  int isDead;
} Vertex;


void Dijkstra(Vertex* graph, int nodecount, int source);

