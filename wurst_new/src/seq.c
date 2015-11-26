/*
 * 26 April 2007
 * Operations on amino sequences of structure type: struct seq
 * rcsid = "$Id: seq.c,v 1.1 2007/06/01 12:22:20 mmosisch Exp $"
 */

#include <string.h>
#include <ctype.h>
#include "e_malloc.h"
#include "read_seq_i.h"
#include "seq.h"
#include "seq_i.h"
#include "str.h"

/* --------------------- seq_duplicate --------------
 * creates a copy of a n-times duplicated sequence.
 * n = 1 returns a sequence copy
 * n = 2 returns a sequence duplicated (simply doubled)
 * n = ...
 * n = 0 is not allowed, returns NULL in failure case
 * params:
 * 		*src - struct seq : sequence to be duplicated
 * 		 n   - size_t     : amount of wanted duplications
 * return:
 * 		dst - struct seq  : new seq
 */
struct seq *
seq_duplicate ( const struct seq *seq, size_t n )
{
	size_t newsize;
    struct seq *dst;

	if (!seq || n < 1)
		return NULL;

	newsize = seq->length * n ;
    dst = E_MALLOC (sizeof (*dst));
    seq_ini (dst);
    dst->seq = E_MALLOC( newsize * sizeof(char));

    if (seq->seq){
    	size_t i;
    	for(i = 0; i < n ; i++){
    		memcpy ( dst->seq + (i * seq->length) , seq->seq , seq->length );
    	}
        /* dst->seq = save_anything ( seq->seq, ( n + 1) * sizeof (seq->seq[0])); */
    }

    if (seq->comment)
        dst->comment = save_str (seq->comment);
    dst->length = newsize;
    dst->format = seq->format;

    return dst;
}
