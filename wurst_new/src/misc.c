/*
 * 25 February 2002
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "misc.h"
#if !defined(lint) && !defined(DONT_SEE_RCS)
  static const char* rcsid = 
  "$Id: misc.c,v 1.2 2005/01/27 11:30:39 torda Exp $";
#endif

/* ---------------- get_nline  --------------------------------
 * Read lines with fgets(), but throw away anything after a
 * comment (#) char.
 * This came from Thomas' get_line(), but it takes another
 * parameter, the size of the input buffer, to prevent overflows.
 * It also knows (and increments) the input line number so the
 * caller can print out errors with line numbers.
 */

char *
get_nline (FILE *fp, char *lbuf, int *nr_line, const size_t maxbuf)
{
    char *pend, *spoint, *s;
    int len;
    while ((s = fgets (lbuf, (int) maxbuf, fp)) != NULL) {
        *nr_line += 1;
        while ( isspace((int)*s ))
            s++;                                /* First non-blank char */
        if ( (spoint = strchr (s, '#')))   /* Find the # comment marker */
            *spoint = '\0';
        if ((len = strlen (s)) == 0)
            continue;
        if ( s [len - 1] == '\n')                 /* Get rid of newline */
            s [--len] = '\0';
        if ( s [0] == '\0')                     /* Hop over blank lines */
            continue;
        for (pend = s + len - 1; pend >= s; pend--)
            if ( isspace ( (int)*pend ) )      /* Delete trailing space */
                *pend = '\0';
            else
                break;
        if (len != 0)
            return s;
    }
    return NULL;
}
