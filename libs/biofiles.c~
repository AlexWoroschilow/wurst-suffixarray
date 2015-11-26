
/*
 *  biofiles.c
 *  helper functions to handle file types
 *  used in bioinformatics
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 07/10/2007 01:56:15 PM CEST
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "stringutils.h"
#include "basic-types.h"
#include "mathematics.h"
#include "biofiles.h"
#include "fileio.h"
#include "charsequence.h"
#include "assert.h"

/*-------------------------------- readfasta ---------------------------------
 *    
 * @brief reads a fasta file 
 * @returns a stringset containing seqs in fastafile
 * @author Steve Hoffmann 
 *   
 */


fasta_t* 
initfasta(void *space) {
  fasta_t *f;

  f = ALLOCMEMORY(space, NULL, fasta_t, 1);
  f->seqs = NULL;
  f->noofseqs = 0;

  return f;
}


void
addfasta(void *space,
    fasta_t *f,
    char    *desc,
    Uint    descrlen,
    char    *sequence,
    Uint    seqlen) {


  f->seqs = ALLOCMEMORY(space, (f->seqs), CharSequence*, (f->noofseqs)+1);
  f->seqs[f->noofseqs] = ALLOCMEMORY(space, NULL, CharSequence, 1);

  f->seqs[f->noofseqs]->description = desc;
  f->seqs[f->noofseqs]->descrlen = descrlen;
  f->seqs[f->noofseqs]->sequence = sequence;
  f->seqs[f->noofseqs]->length = seqlen;

  f->noofseqs++;
  return;
}

fasta_t**
chopfasta(void *space, fasta_t* f, Uint pieces) {
  Uint size, r, i, j, offset=0;
  fasta_t **chops; 

  size = f->noofseqs/pieces;
  r = f->noofseqs%size;

  assert((pieces*size)+r == f->noofseqs);

  chops = ALLOCMEMORY(space, NULL, fasta_t*, pieces);
  for(i=0; i < pieces; i++) {

    chops[i] = ALLOCMEMORY(space, NULL, fasta_t, 1);
    if (i < pieces-1) 
      chops[i]->noofseqs=size;
    else
      chops[i]->noofseqs = size+r;

    chops[i]->seqs = ALLOCMEMORY(space, NULL, CharSequence*, chops[i]->noofseqs);

    for(j=0; j < chops[i]->noofseqs; j++) {  
      chops[i]->seqs[j] = f->seqs[j+offset];
    }

    offset += chops[i]->noofseqs;
  }
  return chops;
}

fasta_t*
readsolexa(void *space, char* filename) {
  stringset_t **s;
  fasta_t *fasta;
  char *descr,
       *sequence;
  Uint i,
       lines;

  s = readcsv(space, filename, "\t", &lines);
  fasta=initfasta(space);

  for(i=0; i < lines; i++) {

    if (s[i]->noofstrings > 4) {
    descr = ALLOCMEMORY(space, NULL, char, 1);
    
    descr[0] = '>';
    descr = concat(space, descr, s[i]->strings[0].str, 1, s[i]->strings[0].len); 
    descr = concat(space, descr, s[i]->strings[1].str, strlen(descr), s[i]->strings[1].len);
    descr = concat(space, descr, s[i]->strings[2].str, strlen(descr), s[i]->strings[2].len); 
    descr = concat(space, descr, s[i]->strings[3].str, strlen(descr), s[i]->strings[2].len); 
    
    sequence = ALLOCMEMORY(space, NULL, char, s[i]->strings[4].len+1);
    memmove(sequence, s[i]->strings[4].str, s[i]->strings[4].len);
    sequence[s[i]->strings[4].len]='\0';

    addfasta(space, fasta, descr, strlen(descr), sequence, strlen(sequence));
    }
    
    destructStringset(space, s[i]);
  }

  FREEMEMORY(space, s);
  return fasta;
}


fasta_t* 
readfasta(void *space, char* filename) {

  char ch;
  char *buffer;
  char *descrbuffer = NULL;
  Uint  descrlength = 0;
  FILE *fp;
  unsigned char desc=1;
  unsigned char first=0;
  fasta_t *fasta;
  Uint buffersize = MAXBUFFERSIZE;
  Uint len=0;  

  buffer = ALLOCMEMORY(space, NULL, char, buffersize);
  fasta = initfasta(space);

  fp = fopen(filename, "r");

  if (fp == NULL) {
    fprintf(stderr,"Couldn't open file '%s'. Exit forced.\n", filename);
    exit(-1);
  }
  
  while((ch=getc(fp)) != EOF) {	

    if(len == buffersize-1) {
      buffersize = 2*buffersize+1;
      buffer = ALLOCMEMORY(space, buffer, char, buffersize);
    }

    if(ch=='>' && !first) {
        desc=1;
        first = 1;
    }

    if(ch=='>' && len > 0) {
      buffer = ALLOCMEMORY(space, buffer, char, len+1);  
      buffer[len]='\0';
      desc = 1;
      addfasta(space, fasta, descrbuffer, descrlength, buffer, len);

      len = 0;
      descrlength = 0;
      descrbuffer = NULL;
      buffersize = MAXBUFFERSIZE;
      buffer = ALLOCMEMORY(space, NULL, char, buffersize);
    }

    if(!desc && ch =='\n') {
      /*do nothing.*/
    } else {
      if(desc && ch == '\n') { 
        buffer = ALLOCMEMORY(space, buffer, char, len+1);  
        buffer[len]='\0'; 

        descrbuffer = buffer;
        descrlength = len;

        len = 0;
        buffersize = MAXBUFFERSIZE;
        buffer = ALLOCMEMORY(space, NULL, char, buffersize);
        desc = 0;
      } else {
        len++;
        buffer[len-1]=(char)ch;
      }
    }
  }

  buffer = ALLOCMEMORY(space, buffer, char, len+1);  
  buffer[len]='\0'; 
  addfasta(space, fasta, descrbuffer, descrlength, buffer, len);

  fclose(fp);
  return fasta;
}


gff_t*
initGff(void *space) {
  gff_t *g;

  g = ALLOCMEMORY(space, NULL, gff_t, 1);
  g->seqname = NULL;
  g->seqnamelen = 0;
  g->source = NULL;
  g->sourcelen = 0;
  g->feat = NULL;
  g->featlen = 0;
  g->start = 0;
  g->end = 0;
  g->score = .0;
  g->strand = '0';
  g->frame = '0';
  g->attrib = NULL;
  g->attriblen = 0;

  return g;
}


void
writeGff(char *filename, gff_t *set, Uint len) {
  Uint i;
  gff_t *g;
  FILE *fp;


  fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "couldn't open %s - exit forced", filename);
    exit(-1);
  }

  for(i=0; i < len; i++) {
    g = &set[i];  
    fprintf(fp,"%s\t%s\t%s\t", g->seqname, g->source, g->feat);
    fprintf(fp,"%d\t%d\t%c\t", g->start, g->end, g->strand);
    fprintf(fp,"%c\t%s\n", g->frame, g->attrib);
  }

  fclose(fp);
  return;
}



