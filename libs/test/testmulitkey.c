
/*
 *  testmulitkey.c
 *  functions to test bentley sedgewick
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/13/06 21:44:33 CET
 *  
 */

 
/*----------------------------------- main -----------------------------------
 *    
 * the main function
 * 
 */
 
int
main (int argc, char** argv)
{
  const char x[]="asimpletreewillbegeneratedfromthisstring$";
  quickSortMultiKey(NULL, x, strlen(x), comparemultikeystring, strlen(x)-1, NULL); 
  
  return EXIT_SUCCESS;
}


