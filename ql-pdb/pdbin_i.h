/*
 * 23 March 2002
 * Interface to pdb reading routines.
 * rcsid = $Id: pdbin_i.h,v 1.1 2002/04/02 00:04:53 torda Exp $
 */
#ifndef PDBIN_I_H
#define PDBIN_I_H

struct coord;
struct coord *
pdb_read ( const char *fname, const char *acq_c, const char chain);

#endif /* PDBIN_I_H */
