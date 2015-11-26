#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "coord.h"
#include "dssp.h"
#include "e_malloc.h"
#include "mprintf.h"
#include "sec_s_i.h"
#include "seq.h"

#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: dssp.c,v 1.13 2005/06/13 11:28:58 torda Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */


enum {SETBITS = 32};


typedef unsigned char boolean;

#ifndef true
enum {
    false = 0,
    true  = 1
};
#endif



/************************************************************
 *                                                           *
 *          p2clib.c                                           *
 *                                                           *
 ************************************************************/

/* Run-time library for use with "p2c", the Pascal to C translator */

/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */

/* File functions */


/* Sets are stored as an array of longs.  S[0] is the size of the set;
   S[N] is the N'th 32-bit chunk of the set.  S[0] equals the maximum
   I such that S[I] is nonzero.  S[0] is zero for an empty set.  Within
   each long, bits are packed from lsb to msb.  The first bit of the
   set is the element with ordinal value 0.  (Thus, for a "set of 5..99",
   the lowest five bits of the first long are unused and always zero.) */

/* (Sets with 32 or fewer elements are normally stored as plain longs.) */

static long *
P_setunion( long *d, long *s1, long *s2)         /* d := s1 + s2 */
{
    long *dbase = d++;
    int sz1 = *s1++, sz2 = *s2++;
    while (sz1 > 0 && sz2 > 0) {
        *d++ = *s1++ | *s2++;
	sz1--, sz2--;
    }
    while (--sz1 >= 0)
	*d++ = *s1++;
    while (--sz2 >= 0)
	*d++ = *s2++;
    *dbase = d - dbase - 1;
    return dbase;
}

static long *
P_setdiff(long *d, long *s1, long *s2)          /* d := s1 - s2 */
{
    long *dbase = d++;
    int sz1 = *s1++, sz2 = *s2++;
    while (--sz1 >= 0 && --sz2 >= 0)
        *d++ = *s1++ & ~*s2++;
    if (sz1 >= 0) {
        while (sz1-- >= 0)
            *d++ = *s1++;
    }
    while (--d > dbase && !*d)
        ;
    *dbase = d - dbase;
    return dbase;
}




static int
P_inset(unsigned val, long *s)                 /* val IN s */
{
    int bit;
    unsigned long *us = (unsigned long *) s;
    bit = val % SETBITS;
    val /= SETBITS;
    if (val < *us++ && ((1<<bit) & us[val]))
	return 1;
    return 0;
}


static long
*P_addset(long *ss, unsigned val)              /* s := s + [val] */
{
    unsigned long *s = (unsigned long *)ss;
    unsigned long *sbase = s;
    unsigned int bit, size;
    bit = val % SETBITS;
    val /= SETBITS;
    size = *s;
    if (++val > size) {
        s += size;
        while (val > size)
            *++s = 0, size++;
        *sbase = size;
    } else
        s += val;
    *s |= 1<<bit;
    return (long *)sbase;
}

static int
P_setequal(long *s1, long *s2)              /* s1 = s2 */
{
    int size = *s1++;
    if (*s2++ != size)
        return 0;
    while (--size >= 0) {
        if (*s1++ != *s2++)
            return 0;
    }
    return 1;
}



static long *
P_setcpy(long *d, long *s)                /* d := s */
{
    long *save_d = d;

#ifdef SETCPY_MEMCPY
    memcpy(d, s, (*s + 1) * sizeof(long));
#else
    int i = *s + 1;
    while (--i >= 0)
        *d++ = *s++;
#endif
    return save_d;
}


/* s is a "smallset", i.e., a 32-bit or less set stored
   directly in a long. */

static long *
P_expset(long *d, long s)                /* d := s */
{
    if (s) {
	d[1] = s;
	*d = 1;
    } else
        *d = 0;
    return d;
}

/* End. */

/*******************  MATHEMATICAL CONSTANTS  **************************
            YVERTEX, - ARE Y,Z-COMPONENTS OF THE FIRST ICOSAHEDRON VERTEX. THE
            ZVERTEX    X-COMPONENT IS 0.
            EPS      - NUMERICAL TOLERANCE
--------------------------------------------------------------------*/

#ifndef M_PI
    static const float M_PI = 3.14159265358979323846;
#endif
#ifndef M_PI_2
    static const float M_PI_2 = 1.57079632679489661923;
#endif

static const float TWOPI   = 2 * 3.14159265358979323846;
static const float RADIAN  = 57.29578;

/***************  ARRAY DIMENSIONING CONSTANTS  ***********************
        NMAX     - MAXIMUM NUMBER OF AMINOACID RESIDUES IN ARRAY CHAIN
        MAXATOM  - MAXIMUM NUMBER OF SIDECHAIN ATOMS IN ARRAY SIDECHAIN
        MAXBRIDGE- MAXIMUM NUMBER OF BRIDGES IN ARRAY BRIDGETABLE
        NFACE,   - NUMBER OF FACES OF POLYHEDRON. THE COORDINATES OF THE CENTRE
        ORDER      OF EACH TRIANGULAR FACE ARE STORED IN ARRAY P, THE AREA
                  IS STORED IN ARRAY WP IN PROCEDURE FLAGACCESS. NFACE MUST
                  BE OF THE FORM NFACE=20*(4**ORDER), ORDER=0,1,2,...
                  THE ACCURACY OF THE SOLVENT ACCESSIBLE SURFACE OF EACH
                  AMINOACID RESIDUE IS ONE ANGSTROM**2 FOR ORDER=2,NFACE=320.
       MAXPACK  - MAXIMUM NUMBER OF PROTEIN ATOMS WHICH CAN INTRUDE INTO
                  SOLVENT AROUND ANY GIVEN TEST ATOM. THE COORDINATES OF
                  THESE ATOMS ARE STORED IN ARRAY X, THEIR RADII IN ARRAY RX
                  IN PROCEDURE SURFACE.
       MAXHIST  - NUMBER OF SLOTS IN ARRAYS HELIXHIST AND BETAHIST USED FOR
                  LENGTH STATISTICS OF SECONDARY STRUCTURE.
       MAXSS    - MAXIMUM NUMBER OF SSBOND RECORDS ON INPUT FILE. THE
                  DISULFIDE BOND ARE SAVED IN ARRAY SSBONDS.
--------------------------------------------------------------------*/

enum {
    NMAX      = 6000,
    MAXATOM   = 40000,
    MAXBRIDGE = 300,
    MAXSIDEATOMS = 20
};


/***/
/*********************  PHYSICAL CONSTANTS   **************************
         SSDIST   - MAXIMUM ALLOWED DISTANCE OF DISULFIDE BRIDGE
         BREAKDIST- MAXIMUM ALLOWED PEPTIDE BOND LENGTH. IF DISTANCE IS
                    GREATER A POLYPEPTIDE CHAIN INTERRUPTION IS ASSUMED.
                    ALL ATOMS OF A RESIDUE
         CADIST   - MINIMUM DISTANCE BETWEEN ALPHA-CARBON ATOMS SUCH THAT NO
                    BACKBONE HYDROGEN BONDS CAN BE FORMED
                    DIST     - SMALLEST ALLOWED DISTANCE BETWEEN ANY ATOMS
         MAXDIST  - LARGEST ALLOWED DISTANCE BETWEEN SIDECHAIN ATOM AND C-ALPHA
                    WITHIN A RESIDUE
         Q        - COUPLING CONSTANT FOR ELECTROstatic ENERGY
                    Q=-332*0.42*0.2*1000.0
         HBLOW    - LOWEST ALLOWED  ENERGY OF A HYDROGEN BOND IN CAL/MOL
         HBHIGH   - HIGHEST ALLOWED ENERGY OF A HYDROGEN BOND IN CAL/MOL
                        --------------------------------------------------------------------*/

static const float SSDIST    = 3.0;
static const float CADIST    = 9.0;
static const float DIST      = 0.5;

static const float Q         = -27888.0;

enum {
    HBLOW = -9900,
    HBHIGH = -500
};

/***/
/***************** GLOBAL DATA TYPE DEFINITIONS ***********************/

typedef float vector[3];
typedef char char4[4];
typedef char char6[6];
typedef enum {
    parallel, antiparallel, nobridge
} bridgetyp;
typedef enum {
    symbol, turn3, turn4, turn5, bend, chirality, beta1, beta2
} structure;

typedef struct hydrogenbond {
    long residue, energy;
} hydrogenbond;

typedef hydrogenbond bonds[2];

typedef struct backbone {
    char6 aaident;
    char sheetlabel, aa;
    char4 threelettercode;
    char ss[(long)beta2 - (long)symbol + 1];
    long partner[(long)beta2 - (long)beta1 + 1];
    /*  long access; */
    float alpha, kappa;
    bonds acceptor, donor;
    vector h, n, ca, c, o;
    long atompointer, nsideatoms;
} backbone;

typedef struct bridge {
    char sheetname, laddername;
    bridgetyp btyp;
    long linkset[MAXBRIDGE / 32 + 2];
    long ib, ie, jb, je, from, towards;
} bridge;

static backbone chain[NMAX + 1];
static vector sidechain[MAXATOM];
static bridge bridgetable[MAXBRIDGE];

static void
VecCopy(float *dest, float *source)
{
    dest[0]=source[0];
    dest[1]=source[1];
    dest[2]=source[2];
}


static float
Atan2(float y, float x)
{
    float z;

    if (x != 0.0)
        z = atan(y / x);
    else if (y > 0.0)
        z = M_PI_2;
    else if (y < 0.0)
        z = -M_PI_2;
    else
        z = TWOPI;
    if (x >= 0.0)
        return z;
    if (y > 0.0)
        z += M_PI;
    else
        z -= M_PI;
    return z;
}  /* Atan2 */


/***/

static void
Diff(float *x, float *y, float *z)
{
    z[0] = x[0] - y[0];
    z[1] = x[1] - y[1];
    z[2] = x[2] - y[2];
}


/***/

static float
Dot(float *x, float *y)
{
    return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}


/***/

static void
Cross(float *x, float *y, float *z)
{
    z[0] = x[1] * y[2] - y[1] * x[2];
    z[1] = x[2] * y[0] - y[2] * x[0];
    z[2] = x[0] * y[1] - y[0] * x[1];
}


/***/



/***/

static float
Dihedralangle(float *v1, float *v2, float *v3, float *v4)
{
    /*CALCULATES TORSION ANGLE OF A SET OF 4 ATOMS V1-V2-V3-V4.
      DIHEDRALANGLE IS THE ANGLE BETWEEN THE PROJECTION OF
      V1-V2 AND THE PROJECTION OF V4-V3 ONTO A PLANE NORMAL TO
      BOND V2-V3.*/
    /***/
    float Result, u, v;
    vector v12, v43, x, y, z, p;

    Diff(v1, v2, v12);
    Diff(v4, v3, v43);
    Diff(v2, v3, z);
    Cross(z, v12, p);
    Cross(z, v43, x);
    Cross(z, x, y);
    u = Dot(x, x);
    v = Dot(y, y);
    Result = 360.0;
    if (u <= 0.0 || v <= 0.0)
        return Result;
    u = Dot(p, x) / sqrt(u);
    v = Dot(p, y) / sqrt(v);
    if (u != 0.0 || v != 0.0)
        return (Atan2(v, u) * RADIAN);
    return Result;
}  /* Dihedralangle */


/***/

static float
Cosangle(float *v1, float *v2, float *v3, float *v4)
{
    vector u, v;
    float x;

    Diff(v1, v2, u);
    Diff(v3, v4, v);
    x = Dot(u, u) * Dot(v, v);
    if (x > 0.0)
        return (Dot(u, v) / sqrt(x));
    else
        return 0.0;
}  /* Cosangle */


/***/

static float
Distance(float *u, float *v)
{
    float temp, temp1, temp2;

    temp = u[0] - v[0];
    temp1 = u[1] - v[1];
    temp2 = u[2] - v[2];
    return sqrt(temp * temp + temp1 * temp1 + temp2 * temp2);
}  /* Distance */


/***/


/*--------------------------------------------------------------------*/

static boolean
Nochainbreak(long i, long j)
{
    long k;
    boolean test;

    test = (i >= 1 && j <= NMAX && i <= j);
    k = i;
    while (test && k <= j) {
        if (chain[k].aa == '!')
            test = false;
        else
            k++;
    }
    return test;
}  /* Nochainbreak */


/***/
/*--------------------------------------------------------------------*/



typedef enum {
    headercard,
    othercard
} cardtype;
/***/


/* Local variables for Inputcoordinates: */
struct LOC_Inputcoordinates {
    long latom, hatoms;
    unsigned *lchain;
    boolean nmissing, camissing, cmissing, omissing, corelimit;
    vector sidecoordinates[MAXSIDEATOMS];
    float dco;
    char4 sideatomnames[MAXSIDEATOMS];
    backbone reszero, resinfo;
} ;

/***/


/* Local variables for Checksideatoms: */
struct LOC_Checksideatoms {
    struct LOC_Inputcoordinates *LINK;
} ;

/***/

/***/

static void
Putresidue(struct LOC_Inputcoordinates *LINK)
{
    /* insert residue into protein chain */
    long i;
    boolean complete;
    long FORLIM;

    complete = !(LINK->nmissing || LINK->camissing || LINK->cmissing ||
                 LINK->omissing);

    LINK->corelimit = (LINK->latom + LINK->resinfo.nsideatoms > MAXATOM ||
                       *LINK->lchain > NMAX - 2);
    if (complete && !LINK->corelimit) {
/*      Checksideatoms(&LINK->resinfo, LINK); */
        VecCopy(LINK->resinfo.h, LINK->resinfo.n);


        if (Nochainbreak(*LINK->lchain, *LINK->lchain) && LINK->resinfo.aa != 'P') {
            LINK->dco = Distance(chain[*LINK->lchain].c, chain[*LINK->lchain].o);
            for (i = 0; i <= 2; i++)
                LINK->resinfo.h[i] = LINK->resinfo.n[i] +
                    (chain[*LINK->lchain].c[i] - chain[*LINK->lchain].o[i]) / LINK->dco;
        }
        (*LINK->lchain)++;
        chain[*LINK->lchain] = LINK->resinfo;
        FORLIM = LINK->resinfo.nsideatoms;
        for (i = 0; i < FORLIM; i++)
            VecCopy(sidechain[LINK->latom + i], LINK->sidecoordinates[i]);
        LINK->latom += LINK->resinfo.nsideatoms;
    }
    if (Nochainbreak(*LINK->lchain, *LINK->lchain) && !complete)
        (*LINK->lchain)++;
    LINK->resinfo = LINK->reszero;
    LINK->nmissing = true;
    LINK->camissing = true;
    LINK->cmissing = true;
    LINK->omissing = true;
}  /* Putresidue */


static void
vcopy ( vector v, struct RPoint *r)
{
    v[0] = r->x;
    v[1] = r->y;
    v[2] = r->z;
}

/***/
/*--------------------------------------------------------------------*/
/* SEE BROOKHAVEN PROTEIN DATA BANK ATOMIC COORDINATE ENTRY FORMAT
   OF DEC. 1981.
   -------------------------------------------------------------------*/

static int
Inputcoordinates(unsigned *lchain_, const struct coord *c)
{
    struct LOC_Inputcoordinates V;
    long i, j;
    size_t k;
    structure s;
/*    cardtype ctype; */
/*  long cardhist[(long)othercard - (long)headercard + 1]; */
    const char *this_sub = "dssp-Inputcoordinates";
    /***/

    V.lchain = lchain_;

    V.latom = 0;
    V.hatoms = 0;
    for (j = 0; j <= 5; j++)
        V.reszero.aaident[j] = ' ';
    V.reszero.aa = '!';
    /*    V.reszero.access = 0; */
    strncpy(V.reszero.threelettercode, "    ", sizeof(char4));
    for (s = symbol; s <= beta2; s++)
        V.reszero.ss[s - symbol] = ' ';
    V.reszero.sheetlabel = ' ';
    V.reszero.partner[0] = 0;
    V.reszero.partner[beta2 - beta1] = 0;
    V.reszero.alpha = 360.0;
    V.reszero.kappa = 360.0;
    for (j = 0; j <= 1; j++) {
        V.reszero.acceptor[j].residue = 0;
        V.reszero.acceptor[j].energy = 0;
        V.reszero.donor[j].residue = 0;
        V.reszero.donor[j].energy = 0;
    }
    V.reszero.atompointer = 0;
    V.reszero.nsideatoms = 0;
    for (j = 0; j <= 2; j++) {
        V.reszero.h[j] = 0.0;
        V.reszero.n[j] = 0.0;
        V.reszero.ca[j] = 0.0;
        V.reszero.c[j] = 0.0;
        V.reszero.o[j] = 0.0;
    }
    for (i = 0; i <= NMAX; i++)
        chain[i] = V.reszero;
/*  for (ctype = headercard; ctype <= othercard;ctype++)
 *      cardhis[ctype - headercard] = 0;
 */
    V.corelimit = false;

    V.resinfo = V.reszero;
    V.nmissing = true;
    V.camissing = true;
    V.cmissing = true;
    V.omissing = true;

    for (k = 0; k < c->size; k++) {
        vcopy (V.resinfo.n,  c->rp_n + k);
        vcopy (V.resinfo.ca, c->rp_ca + k);
        vcopy (V.resinfo.c,  c->rp_c + k);
        vcopy (V.resinfo.o,  c->rp_o + k);
        /* Don't worry about c-beta */
        V.camissing = V.nmissing = V.cmissing = V.omissing = false;
        V.resinfo.aa = toupper (c->seq->seq[k]);
        Putresidue (&V);
    }

    if (V.corelimit) {
        err_printf(this_sub, " !!! NUMBER OF ATOMS OR RESIDUES EXCEEDS ");
        err_printf(this_sub, "STORAGE CAPACITY !!!\n");
    }
    if (!Nochainbreak(*V.lchain, *V.lchain))
        (*V.lchain)--;

    if (*V.lchain < 1) {
        err_printf(this_sub, " !!! NO RESIDUE WITH COMPLETE BACKBONE !!!\n");
        return EXIT_FAILURE;
    }

    return (EXIT_SUCCESS);
}  /* Inputcoordinates */



/***/
/*--------------------------------------------------------------------*/

static boolean
Testbond(long i, long j)
{
    /* TESTBOND IS TRUE IF I IS DONOR[=NH] TO J, OTHERWISE FALSE */
    backbone *WITH;

    WITH = chain + i;
    return ((WITH->acceptor[0].residue == j && WITH->acceptor[0].energy < HBHIGH) ||
            (WITH->acceptor[1].residue == j && WITH->acceptor[1].energy < HBHIGH));

}  /* Testbond */


/***/


/*--------------------------------------------------------------------*/

static void
Flagssbonds(const unsigned lchain)
{
    boolean ssbond;

    long i, ii, jj;
    unsigned j;

    long FORLIM;

    /***/



    FORLIM = lchain - 2;
    for (i = 1; i <= FORLIM; i++) {
        if (chain[i].aa == 'C' && chain[i].nsideatoms > 1) {
            ii = chain[i].atompointer + 2;
            j = i + 1;
            do {
                j++;
                ssbond = false;
                if (chain[j].nsideatoms > 1 && chain[j].aa == 'C')
                    jj = chain[j].atompointer + 2;
                else
                    jj = 0;
                if (jj > 0)
                    ssbond = (Distance(sidechain[ii - 1], sidechain[jj - 1]) < SSDIST);
            } while (!(ssbond || j == lchain));

        }
    }

}  /* Flagssbonds */


/***/
/*--------------------------------------------------------------------*/

static void
Flagchirality( const unsigned lchain)
{
    long i;
    float ckap, skap;
    long FORLIM;
    backbone *WITH;

    FORLIM = lchain - 2;
    for (i = 2; i <= FORLIM; i++) {
        WITH = &chain[i];
        if (Nochainbreak(i - 1, i + 2)) {
            WITH->alpha = Dihedralangle(chain[i - 1].ca, WITH->ca, chain[i + 1].ca,
                                        chain[i + 2].ca);
            if (WITH->alpha < 0.0)
                WITH->ss[chirality - symbol] = '-';
            else
                WITH->ss[chirality - symbol] = '+';
        }
    }
    FORLIM = lchain - 2;
    /***/
    for (i = 3; i <= FORLIM; i++) {
        WITH = &chain[i];
        if (Nochainbreak(i - 2, i + 2)) {
            ckap = Cosangle(chain[i].ca, chain[i - 2].ca, chain[i + 2].ca,
                            chain[i].ca);
            skap = sqrt(1.0 - ckap * ckap);
            WITH->kappa = RADIAN * Atan2(skap, ckap);
        }
    }
}  /* Flagchirality */


/***/

static long
Bondenergy(long i, long j)
{
    /*RESIDUE I IS DONOR[=NH],J IS ACCEPTOR[=CO] OF THE PROTON IN THE
      HYDROGEN BOND. THE BONDENERGY IS IN CAL/MOL */
    float dho, dhc, dnc, dno;
    long hbe;
    backbone *WITH;

    hbe = 0;
    WITH = &chain[i];
    if (WITH->aa == 'P')
        return hbe;
    dho = Distance(WITH->h, chain[j].o);
    dhc = Distance(WITH->h, chain[j].c);
    dnc = Distance(WITH->n, chain[j].c);
    dno = Distance(WITH->n, chain[j].o);
    if (dho < DIST || dhc < DIST || dnc < DIST || dno < DIST)
        hbe = HBLOW;
    else
        hbe = (long)floor(Q / dho - Q / dhc + Q / dnc - Q / dno + 0.5);
    if (hbe > HBLOW)
        return hbe;
    hbe = HBLOW;
    return hbe;
}  /* Bondenergy */

/***/

static void
Updatebonds( hydrogenbond *b, hydrogenbond hb)
{
    if (hb.energy < b[0].energy) {
        b[1] = b[0];
        b[0] = hb;
    } else if (hb.energy < b[1].energy)
        b[1] = hb;
}  /* Updatebonds */

/***/

static void
Setbonds(long i, long j)
{
    /*I IS NH, J IS CO*/
    hydrogenbond hb;

    hb.energy = Bondenergy(i, j);
    hb.residue = j;
    /* CO(J) IS ACCEPTOR OF NH(I) */
    Updatebonds(chain[i].acceptor, hb);
    hb.residue = i;
    Updatebonds(chain[j].donor, hb);
}  /* Setbond */


/***/
/*--------------------------------------------------------------------*/

static void
Flaghydrogenbonds( const unsigned lchain )
{
    long i, j, FORLIM;
    backbone *WITH;
    long FORLIM1;

    /***/

    FORLIM = lchain;
    for (i = 1; i <= FORLIM; i++) {
        if (Nochainbreak(i, i)) {
            WITH = &chain[i];
            FORLIM1 = lchain;
            for (j = i + 1; j <= FORLIM1; j++) {
                if (Nochainbreak(j, j)) {
                    if (Distance(WITH->ca, chain[j].ca) < CADIST) {
                        Setbonds(i, j);
                        if (j != i + 1)
                            Setbonds(j, i);
                    }
                }
            }
        }
    }
}  /* Flaghydrogenbonds */


/***/

static int
Ladder( long i, long j, bridgetyp b, unsigned *nbridge)
{
    long k;
    boolean found;
    bridge *WITH;

    k = 1;
    if (b == nobridge || i >= j)
        return EXIT_SUCCESS;
    do {
        WITH = &bridgetable[k - 1];
        if (WITH->ib == 0) {
            WITH->ib = i;
            WITH->ie = i;
            WITH->jb = j;
            WITH->je = j;
            WITH->from = 0;
            WITH->towards = 0;
            WITH->btyp = b;
            (* nbridge)++;
            found = true;
        } else {
            found = (WITH->btyp == b && i == WITH->ie + 1) & Nochainbreak(WITH->ie,i) &
                (((j == WITH->je + 1 && b == parallel) & Nochainbreak(WITH->je, j)) |
                 ((j == WITH->jb - 1 && b == antiparallel) & Nochainbreak(j, WITH->jb)));

            if (found) {
                WITH->ie++;
                if (b == parallel)
                    WITH->je++;
                else
                    WITH->jb--;
            } else {
                k++;
                if (k > MAXBRIDGE) {
                    err_printf("dssp-ladder", "!!! BRIDGETABLE OVERFLOW !!\n");
                    return EXIT_FAILURE;
                }
            }
        }
    } while (!found);   /* Ladder */
    return EXIT_SUCCESS;
}

/***/

static int
Testbridge(long i, const unsigned lchain, unsigned *nbridge)
{
    unsigned j, j2, jj1;
    bridgetyp b;


    j2 = 0;
    j = i + 3;
    if (!Nochainbreak(i - 1, i + 1))
        return (EXIT_SUCCESS);
    jj1 = 0;
    while (j2 == 0 && j < lchain) {
        if (Nochainbreak(j - 1, j + 1)) {
            if ((Testbond(i + 1, j) & Testbond(j, i - 1)) |
                (Testbond(j + 1, i) & Testbond(i, j - 1)))
                b = parallel;
            else if ((Testbond(i + 1, j - 1) & Testbond(j + 1, i - 1)) |
                     (Testbond(j, i) & Testbond(i, j)))
                b = antiparallel;
            else
                b = nobridge;
            if (b != nobridge) {
                if (jj1 == 0) {
                    jj1 = j;
                    if (Ladder(i, j, b, nbridge) == EXIT_FAILURE)
                        return (EXIT_FAILURE);
                } else if (j != jj1) {
                    j2 = j;
                    if (Ladder(i, j, b, nbridge) == EXIT_FAILURE)
                        return (EXIT_FAILURE);
                }
            }
        }
        j++;
    }
    return (EXIT_SUCCESS);
}  /* Testbridge */

/***/

static void
Extendladder( const unsigned nbridge)
{
    long i, ib1, jb1, je1;
    unsigned j;
    boolean bulge;
    long FORLIM;
    bridge *WITH;
    const char *this_sub = "dssp-Extendladder";

    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
        WITH = &bridgetable[i - 1];
        j = i + 1;
        while (j <= nbridge && WITH->towards == 0) {
            ib1 = bridgetable[j - 1].ib;
            jb1 = bridgetable[j - 1].jb;
            je1 = bridgetable[j - 1].je;
            bulge = (Nochainbreak(WITH->ie, ib1) && ib1 - WITH->ie < 6 &&
                     bridgetable[j - 1].btyp == WITH->btyp &&
                     bridgetable[j - 1].from == 0);
            if (bulge) {
                switch (WITH->btyp) {

                case parallel:
                    bulge = ((jb1 - WITH->je < 6 && ib1 - WITH->ie < 3) ||
                             jb1 - WITH->je < 3) & Nochainbreak(WITH->je, jb1);
                    break;

                case antiparallel:
                    bulge = ((WITH->jb - je1 < 6 && ib1 - WITH->ie < 3) ||
                             WITH->jb - je1 < 3) & Nochainbreak(je1, WITH->jb);
                    break;
                case nobridge:
                    err_printf (this_sub, "Bug %s %d\n", __FILE__, __LINE__);
                    exit (EXIT_FAILURE);
                }
            }
            if (bulge) {
                WITH->towards = j;
                bridgetable[j - 1].from = i;
            }
            j++;
        }
    }
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
        WITH = &bridgetable[i - 1];
        if (WITH->from == 0) {
            P_expset(WITH->linkset, 0L);
            j = i;
            do {
                P_addset(WITH->linkset, (int)j);
                j = bridgetable[j - 1].towards;
            } while (j != 0);
            j = WITH->towards;
            while (j != 0) {
                P_setcpy(bridgetable[j - 1].linkset, WITH->linkset);
                j = bridgetable[j - 1].towards;
            }
        }
    }
}  /* Extendladder */

/* Local variables for Sheet: */
struct LOC_Sheet {
    long ladderset[MAXBRIDGE / 32 + 2], sheetset[MAXBRIDGE / 32 + 2];
} ;

/***/

static boolean
Link(long l1, long l2)
{
    /* LINK IS TRUE IF THERE IS A COMMON RESIDUE IN LADDERS L1 AND L2 */
    long ib1, ie1, jb1, je1, ib2, ie2, jb2, je2;

    ib1 = bridgetable[l1 - 1].ib;
    ie1 = bridgetable[l1 - 1].ie;
    jb1 = bridgetable[l1 - 1].jb;
    je1 = bridgetable[l1 - 1].je;
    ib2 = bridgetable[l2 - 1].ib;
    ie2 = bridgetable[l2 - 1].ie;
    jb2 = bridgetable[l2 - 1].jb;
    je2 = bridgetable[l2 - 1].je;
    return ((ie1 >= ib2 && ib1 <= ie2) || (ie1 >= jb2 && ib1 <= je2) ||
            (je1 >= ib2 && jb1 <= ie2) || (je1 >= jb2 && jb1 <= je2));
}  /* Link */

/***/

static void
Findsheet( struct LOC_Sheet *LINK, const unsigned nbridge)
{
    unsigned l1;
    boolean finish;
    unsigned FORLIM, FORLIM1;

    /***/

    P_expset(LINK->sheetset, 0L);
    l1 = 0;
    if (*LINK->ladderset != 0L) {
        do {
            l1++;
        } while (!P_inset((int)l1, LINK->ladderset));
    }
    if (l1 > 0)
        P_setcpy(LINK->sheetset, bridgetable[l1 - 1].linkset);
    if (l1 == 0)
        return;
    do {
        unsigned l2, l3;
        finish = true;
        FORLIM = nbridge;
        for (l3 = 1; l3 <= FORLIM; l3++) {
            if (P_inset((int)l3, LINK->sheetset)) {
                FORLIM1 = nbridge;
                for (l2 = 1; l2 <= FORLIM1; l2++) {
                    if (P_inset((int)l2, LINK->ladderset)) {
                        if (Link(l3, l2)) {
                            P_setunion(LINK->sheetset, LINK->sheetset,
                                       bridgetable[l2 - 1].linkset);
                            P_setdiff(LINK->ladderset, LINK->ladderset,
                                      bridgetable[l2 - 1].linkset);
                            finish = false;
                        }
                    }
                }
            }
        }
    } while (!finish);   /* Findsheet */
}

/***/

static void
Sheet(const unsigned nbridge )
{
    struct LOC_Sheet V;
    long asci, i, j;
    char ccs;
    long FORLIM;
    bridge *WITH;

    P_expset(V.ladderset, 0L);
    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++)
        P_addset(V.ladderset, (int)i);
    ccs = '@';
    asci = 64;
    while (*V.ladderset != 0L) {
        ccs++;
        if (ccs > 'z') {
            ccs = 'A';
        }
        Findsheet(&V, nbridge);
        FORLIM = nbridge;
        for (i = 1; i <= FORLIM; i++) {
            WITH = &bridgetable[i - 1];
            if (P_inset((int)i, V.sheetset) && WITH->from == 0) {
                if (asci == 90) {
                    asci = 64;
                }
                asci++;
                if (WITH->btyp == parallel)
                    WITH->laddername = (char)(asci + 32);
                else
                    WITH->laddername = (char)asci;
                WITH->sheetname = ccs;
                P_setcpy(WITH->linkset, V.sheetset);
                j = WITH->towards;
                while (j != 0) {
                    bridgetable[j - 1].laddername = WITH->laddername;
                    bridgetable[j - 1].sheetname = WITH->sheetname;
                    P_setcpy(bridgetable[j - 1].linkset, V.sheetset);
                    j = bridgetable[j - 1].towards;
                }
            }
        }
    }
}  /* Sheet */

/***/

static void
Markstrands( const unsigned lchain, const unsigned nbridge)
{
    long i, j, l, ib0, ie0, jb0, je0;
    structure beta, betai, betaj;
    long iset[beta2 - beta1 + 1][9],
        jset[beta2 - beta1 + 1][9];
    char cc;
    long FORLIM, FORLIM1;
    long SET[9];
    bridge *WITH;
    backbone *WITH1;
    long SET1[9];
    long SET2[3];


    FORLIM = nbridge;
    for (i = 1; i <= FORLIM; i++) {
        if (bridgetable[i - 1].from == 0) {
            j = i;
            for (beta = beta1; beta <= beta2; beta++) {
                P_setcpy(iset[beta - beta1], P_expset(SET, 0L));
                P_setcpy(jset[beta - beta1], P_expset(SET, 0L));
            }
            ib0 = lchain;
            ie0 = 0;
            jb0 = lchain;
            je0 = 0;
            do {
                WITH = &bridgetable[j - 1];
                FORLIM1 = WITH->ie;
                for (l = WITH->ib; l <= FORLIM1; l++) {
                    WITH1 = &chain[l];
                    for (beta = beta1; beta <= beta2; beta++)
                        P_setcpy(iset[beta - beta1], P_setunion(SET1,
                                                                            iset[beta - beta1],
                                                                            P_addset(P_expset(SET, 0L),
                                                                                     WITH1->ss[beta - symbol])));
                }
                FORLIM1 = WITH->je;
                for (l = WITH->jb; l <= FORLIM1; l++) {
                    WITH1 = &chain[l];
                    for (beta = beta1; beta <= beta2; beta++)
                        P_setcpy(jset[beta - beta1], P_setunion(SET1,
                                                                jset[beta - beta1],
                                                                P_addset(P_expset(SET, 0L),
                                                                                     WITH1->ss[beta - symbol])));
                }
                if (WITH->ib < ib0)
                    ib0 = WITH->ib;
                if (WITH->ie > ie0)
                    ie0 = WITH->ie;
                if (WITH->jb < jb0)
                    jb0 = WITH->jb;
                if (WITH->je > je0)
                    je0 = WITH->je;
                j = WITH->towards;
            } while (j != 0);
            j = i;
            if (P_setequal(iset[0], P_addset(P_expset(SET2, 0L), ' ')))
                betai = beta1;
            else
                betai = beta2;
            if (P_setequal(jset[0], P_addset(P_expset(SET2, 0L), ' ')))
                betaj = beta1;
            else
                betaj = beta2;
            do {
                WITH = &bridgetable[j - 1];
                FORLIM1 = WITH->ie;
                for (l = WITH->ib; l <= FORLIM1; l++) {
                    WITH1 = &chain[l];
                    WITH1->ss[betai - symbol] = WITH->laddername;
                    if (WITH->btyp == parallel)
                        WITH1->partner[betai - beta1] = WITH->jb + l - WITH->ib;
                    else
                        WITH1->partner[betai - beta1] = WITH->je - l + WITH->ib;
                }
                FORLIM1 = WITH->je;
                for (l = WITH->jb; l <= FORLIM1; l++) {
                    WITH1 = &chain[l];
                    WITH1->ss[betaj - symbol] = WITH->laddername;
                    if (WITH->btyp == parallel)
                        WITH1->partner[betaj - beta1] = WITH->ib + l - WITH->jb;
                    else
                        WITH1->partner[betaj - beta1] = WITH->ie - l + WITH->jb;
                }
                j = WITH->towards;
            } while (j != 0);
            if (ib0 == ie0)
                cc = 'B';
            else
                cc = 'E';
            for (j = ib0; j <= ie0; j++) {
                WITH1 = &chain[j];
                if (WITH1->ss[0] != 'E')
                    WITH1->ss[0] = cc;
            }
            for (j = jb0; j <= je0; j++) {
                WITH1 = &chain[j];
                if (WITH1->ss[0] != 'E')
                    WITH1->ss[0] = cc;
            }
        }
    }
    FORLIM = nbridge;
    for (j = 0; j < FORLIM; j++) {
        WITH = &bridgetable[j];
        FORLIM1 = WITH->ie;
        for (l = WITH->ib; l <= FORLIM1; l++)
            chain[l].sheetlabel = WITH->sheetname;
        FORLIM1 = WITH->je;
        for (l = WITH->jb; l <= FORLIM1; l++)
            chain[l].sheetlabel = WITH->sheetname;
    }
}  /* Markstrands */


/***/
/*--------------------------------------------------------------------*/

static int
Flagbridge( const unsigned lchain )
{
    long i, FORLIM;
    bridge *WITH;
    unsigned nbridge = 0;

    for (i = 0; i < MAXBRIDGE; i++) {
        WITH = &bridgetable[i];
        WITH->ib = 0;
        WITH->ie = 0;
        WITH->jb = 0;
        WITH->je = 0;
        WITH->btyp = nobridge;
    }

    FORLIM = lchain;
    for (i = 2; i < FORLIM; i++)
        if (Testbridge(i, lchain, &nbridge) == EXIT_FAILURE)
            return (EXIT_FAILURE);
    if (nbridge == 0)
        return (EXIT_SUCCESS);
    Extendladder(nbridge);
    Sheet(nbridge);
    Markstrands(lchain, nbridge);
    return (EXIT_SUCCESS);
}  /* Flagbridge */


/*--------------------------------------------------------------------*/
static void
Flagsymbol(const unsigned lchain)
{
    /* FLAGS ALPHA HELICES AND TURNS IN SYMBOL COLUMN */
    long i, j, k;
    char cc;
    long nhset[9];
    structure turn;
    boolean empty;
    long FORLIM;
    backbone *WITH;

    P_addset(P_expset(nhset, 0L), '>');
    P_addset(nhset, 'X');
    FORLIM = lchain - 4;
    for (i = 2; i <= FORLIM; i++) {
        if (P_inset(chain[i - 1].ss[turn4 - symbol], nhset) &
            P_inset(chain[i].ss[turn4 - symbol], nhset)) {
            for (j = i; j <= i + 3; j++)
                chain[j].ss[0] = 'H';
        }
    }
    FORLIM = lchain - 3;
    for (i = 2; i <= FORLIM; i++) {
        if (P_inset(chain[i - 1].ss[turn3 - symbol], nhset) &
            P_inset(chain[i].ss[turn3 - symbol], nhset)) {
            empty = true;
            for (j = i; j <= i + 2; j++) {
                WITH = &chain[j];
                if (WITH->ss[0] != 'G' && WITH->ss[0] != ' ')
                    empty = false;
            }
            if (empty) {
                for (j = i; j <= i + 2; j++)
                    chain[j].ss[0] = 'G';
            }
        }
    }
    FORLIM = lchain - 5;
    for (i = 2; i <= FORLIM; i++) {
        if (P_inset(chain[i - 1].ss[turn5 - symbol], nhset) &
            P_inset(chain[i].ss[turn5 - symbol], nhset)) {
            empty = true;
            for (j = i; j <= i + 4; j++) {
                WITH = &chain[j];
                if (WITH->ss[0] != 'I' && WITH->ss[0] != ' ')
                    empty = false;
            }
            if (empty) {
                for (j = i; j <= i + 4; j++)
                    chain[j].ss[0] = 'I';
            }
        }
    }
    FORLIM = lchain;
    for (i = 2; i < FORLIM; i++) {
        WITH = &chain[i];
        if (WITH->ss[0] == ' ') {
            cc = ' ';
            j = 1;
            for (turn = turn3; turn <= turn5; turn++) {
                j++;
                for (k = 1; k <= j; k++) {
                    if (i > k) {
                        if (P_inset(chain[i - k].ss[turn - symbol], nhset))
                            cc = 'T';
                    }
                }
            }
            if (cc == ' ')
                cc = WITH->ss[bend - symbol];
            WITH->ss[0] = cc;
        }
    }
}  /* Flagsymbol */

/* ---------------- Flagturn ----------------------------------
 */
static void
Flagturn( const unsigned lchain )
{
    long i, j, k;
    structure turn;
    char cc;
    long FORLIM1;
    backbone *WITH;

    /***/

    k = 2;
    cc = '2';
    for (turn = turn3; turn <= turn5; turn++) {
        k++;
        cc++;
        FORLIM1 = lchain - k;
        for (i = 1; i <= FORLIM1; i++) {
            if (Nochainbreak(i, i + k)) {
                if (Testbond(i + k, i)) {
                    chain[i + k].ss[turn - symbol] = '<';
                    for (j = 1; j < k; j++) {
                        WITH = &chain[i + j];
                        if (WITH->ss[turn - symbol] == ' ')
                            WITH->ss[turn - symbol] = cc;
                    }
                    WITH = &chain[i];
                    if (WITH->ss[turn - symbol] == '<')
                        WITH->ss[turn - symbol] = 'X';
                    else
                        WITH->ss[turn - symbol] = '>';
                }
            }
        }
    }
    FORLIM1 = lchain;
    for (i = 1; i <= FORLIM1; i++) {
        WITH = &chain[i];
        if (WITH->kappa != 360.0 && WITH->kappa > 70.0)
            WITH->ss[bend - symbol] = 'S';
    }
    Flagsymbol(lchain);
}  /* Flagturn */

/* ---------------- Printout ----------------------------------
 */
static int
Printout(const unsigned lchain, struct coord *c )
{
    backbone *b;
    size_t k;
    const char *this_sub = "dssp-Printout";
    if (lchain != c->size) {
        err_printf (this_sub, "lchain != c->size, %d != %ld\n",
                    lchain, (long int) c->size);
        return EXIT_FAILURE;
    }
    b = chain + 1;
    if (c->sec_typ == NULL)
        c->sec_typ = E_MALLOC (c->size * sizeof (c->sec_typ[0]));
    for (k = 0; k < c->size; k++, b++)
        c->sec_typ[k] = char2ss(b->ss[0]);
#   ifdef Debug_till_i_fall_over_and_bleed
        for (i = 1; i <= lchain; i++) {
            WITH = chain + i;
            mprintf("%5ld ", i);
            for (j = 0; j <= 5; j++)
                mputc(WITH->aaident[j], stdout);
            mprintf(" %c  %c ", WITH->aa, WITH->ss[0]);
            mprintf ("\n");
        }
#   endif
    return (EXIT_SUCCESS);
}

/* ---------------- dssp    -----------------------------------
 */

int
dssp ( struct coord *c)
{
    unsigned lchain = 0;
    if (Inputcoordinates(&lchain, c) == EXIT_FAILURE)
        return EXIT_FAILURE;

    Flagssbonds(lchain);

    Flagchirality(lchain);

    Flaghydrogenbonds(lchain);

    if (Flagbridge(lchain) == EXIT_FAILURE)
        return (EXIT_FAILURE);

    Flagturn(lchain);

    if (Printout(lchain, c) == EXIT_FAILURE)
        return (EXIT_FAILURE);

    return EXIT_SUCCESS;

}
