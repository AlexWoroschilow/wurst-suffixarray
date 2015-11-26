/*
 * 30 May 2002
 * We need a version marker on our coordinate files.
 * We want a revision number, but if we include the string in our
 * file, coord.c, the version control system will grab hold of
 * it.
 */

#include "binver.h"

/* ---------------- bin_version   -----------------------------
 */
const char *
bin_version ( void )
{
    return "wurst $Revision: 1.9 $ version\n";
}
