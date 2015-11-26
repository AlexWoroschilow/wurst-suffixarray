#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "dssp.h"

typedef struct coord Coord;
MODULE = Dssp		PACKAGE = Dssp		
PROTOTYPES: ENABLE

int
dssp(c)
        Coord *c
    CODE:
        RETVAL = dssp (c);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;        
