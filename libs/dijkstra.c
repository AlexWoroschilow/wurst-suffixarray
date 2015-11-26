
/*
 *  dijkstra.c
 *  the famous algorithm
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 09/15/2007 03:29:48 PM CEST
 *  
 */


void Dijkstra(Vertex* graph, int nodecount, int source) {
  for(int i = 0; i < nodecount; i++) {
    if(i == source) {
      graph[i].distance = 0;
      graph[i].isDead = 0;
    } else {
      graph[i].distance = INFINITY;
      graph[i].isDead = 0;
    }
  }
  for(int i = 0; i < nodecount; i++) {
    int next;
    int min = INFINITY+1;
    for(int j = 0; j < nodecount; j++) {
      if(!graph[j].isDead && graph[j].distance < min) {
        next = j;
        min = graph[j].distance;
      }
    }
    for(int j = 0; j < graph[next].numconnect; j++) {
      if(graph[graph[next].connections[j].dest].distance >
          graph[next].distance + graph[next].connections[j].weight)
      {
        graph[graph[next].connections[j].dest].distance =
          graph[next].distance + graph[next].connections[j].weight;
      }
    }
    graph[next].isDead = 1;
  }
  for(int i = 0; i < nodecount; i++) {
    printf("The distance between nodes %i and %i is %i",
        source, i, graph[i].distance);
  }
}
