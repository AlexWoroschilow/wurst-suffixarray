#ifndef BITARRAY_H
#define BITARRAY_H

/*
 *
 *	bitArray.h
 *  declarations for bit arrays
 * 
 *  @author Steve Hoffmann, steve@bioinf.uni-leipzig.de
 *  @company Bioinformatics, University of Leipzig 
 *  @date 07/14/2007 04:15:27 PM CEST  
 *
 */

typedef unsigned char* bitarray;

bitarray initbitarray(void *, Uint length);
void dumpbitarray(bitarray a, Uint len);
unsigned char valbitarray(bitarray a, Uint len, unsigned char val);
void setbitarray(bitarray a, Uint len, unsigned char val);
bitarray resizebitarray(void *space, bitarray, Uint len);
static inline void
setbit(bitarray a, Uint pos, unsigned char val) {
  int  bytes,
       bits,
       shift;
  unsigned char byte, resc;

  bytes = pos >> 3;
  bits  = pos & 7;
  resc = a[bytes];
  
  byte = (unsigned char) val;
  shift = 7 ^ bits;

  byte = byte << shift;
  byte = byte ^ resc;
  byte = byte & (1 << shift);
  
  a[bytes] = byte ^ resc;
}

static inline unsigned char
getbit(bitarray a, Uint pos) {
  int bytes,
      bits;
  unsigned char byte;

  bytes = pos >> 3;
  bits  = pos & 7;
    
  byte = a[bytes];
  byte = byte >> (7 ^ bits);
  byte = byte << 7;
  byte = byte >> 7;
  byte = byte & 1;

  return byte;
}


#endif
