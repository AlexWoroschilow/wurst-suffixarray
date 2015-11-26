#ifndef SUF_ARRAY_H
#define SUF_ARRAY_H

/*
 *
 *	sufarray.h
 *  declarations for enhanced suffix arrays
 *  for char alphabets
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/10/06 22:01:33 CET  
 *
 */

#include "basic-types.h"
#include "falphabet.h"
#include "charsequence.h"
#include "list.h"
#include "multicharseq.h"
#include "openssl/md5.h"

#define LCP_TAB_STORED      ((unsigned char) (1 << 0))
#define CHLD_TAB_STORED     ((unsigned char) (1 << 1))
#define SUFLINK_TAB_STORED  ((unsigned char) (1 << 2))
#define SUFLINK_COMPRESSED  ((unsigned char) (1 << 3))


typedef struct {
    Uint *suf;
    Uint pos;
} suffix_t;

typedef struct {
 /* int up;
    int down;
    int nextlIndex;
 */

  int val;
} childtab;


typedef struct {
  MultiCharSeq  *seq;

  Uint		numofsuffixes;
  char		**suffixptr;
  Uint 		*suftab;
  Uint		*inv_suftab;
  Uint      *suflink;
  Uint      *suflink_l;
  Uint      *suflink_r;

  Uint 		*bwttab; 	     /* burrows-wheeler array*/
  Uint		*lcptab;         /* alternative: Abouelhoda et al.*/

  unsigned char *lcpctab;    /* nB to store lcp values < 255*/
  PairUint	    *llvtab;     /* array of 8B to store lcp val >=255*/
  Uint           llvcnt;
  Uint           maxlcp;

  signed char   *id;
  PairSint      *idvtab;
  Uint          idvcnt;

  childtab      *chldtab;    /* a child table*/
  Uint		    *bcktab;     /* the bucket container*/

} Suffixarray;


Suffixarray* readSuffixarray(void *, char *, CharSequence **, Uint);
void writeSuffixarray(Suffixarray *s, char *filename);
Suffixarray* constructSufArr(void *, CharSequence **, Uint, FAlphabet *);
void constructchildtab(void *, Suffixarray *);
void constructsuflinks(void *, Suffixarray *, Uint *);
void constructLcp (void *, Suffixarray *);
void computeId(void*, Suffixarray *);
Uint* getsufsucc(void *, Suffixarray *);

void destructSufArr (void *, Suffixarray *); 
List* getChildintervals(void *, Suffixarray *, Uint, Uint); 
PairUint getCharInterval(void *, Suffixarray *, Uint, Uint, Uint, char);

inline PairUint getSuflink(Suffixarray *, Uint, Uint);
inline Uint getlcpval(Suffixarray *, Uint, Uint);
inline Uint getfirstlindex(Suffixarray *, Uint, Uint);
inline int id (Suffixarray *, Uint);
void destructinterval(void *space, void *data);

void dumplcps (Suffixarray *);

#endif

