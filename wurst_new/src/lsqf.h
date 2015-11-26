/* "$Id: lsqf.h,v 1.1 2006/02/23 14:03:45 tmargraf Exp $" */
#ifndef LSQF_H 
#define LSQF_H

int
coord_rmsd (struct pair_set *pairset, struct coord *coord1, struct coord *coord2,
	    int sub_flag, float *rmsd, struct coord **c1_new, struct coord **c2_new);

int get_rmsd(struct pair_set *pairset, struct coord *r1,
       struct coord *r2, float *rmsd_ptr, int *count);
#endif
