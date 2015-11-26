#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "align_i.h"
#include "altscores.h"
#include "bad_angle.h"
#include "class_model.h"
#include "classifyStructure.h"
#include "cmp_dmat_i.h"
#include "coord.h"
#include "coord_i.h"
#include "coordinfo_i.h"
#include "dipep_i.h"
#include "geo_gap.h"
#include "lsqf.h"
#include "model.h"
#include "multialign.h"
#include "pair_set_chim.h"
#include "pair_set_i.h"
#include "pair_set_p_i.h"
#include "param_fx_i.h"
#include "pdbin_i.h"
#include "pdbout_i.h"
#include "prob_vec_i.h"
#include "read_ac.h"
#include "read_ac_i.h"
#include "read_ac_strct.h"
#include "read_ac_strct_i.h"
#include "read_blst.h"
#include "read_mat.h"
#include "read_sec_i.h"
#include "read_seq_i.h"
#include "rescore.h"
#include "score_fx_i.h"
#include "score_mat_i.h"
#include "score_probvec.h"
#include "score_sec_i.h"
#include "score_smat.h"
#include "scratch.h"
#include "sec_s_i.h"
#include "silly.h"
#include "scor_set.h"
#include "scor_set_i.h"
#include "seq.h"
#include "seq_i.h"

static int
not_here(char *s)
{
    croak("%s not implemented on this architecture", s);
    return -1;
}

static double
constant(char *name, int len, int arg)
{
    errno = 0;
    if (strEQ(name + 0, "WURST_H")) {    /*  removed */
#       ifdef WURST_H
        return WURST_H;
#       else
        errno = ENOENT;
        return 0;
#       endif
    }
    errno = EINVAL;
    return 0;
}

typedef struct aa_clssfcn       Aa_clssfcn;
typedef struct aa_strct_clssfcn Aa_strct_clssfcn ;
typedef struct dpt_list         Dpt_list;
typedef struct prob_vec         Prob_vec;
typedef struct sec_s_data       Sec_s_data;
typedef struct seq              Seq;
typedef struct seq_array        Seq_array;
typedef struct sub_mat          Sub_mat;
typedef struct score_mat        Score_mat;
typedef struct pair_set         Pair_set;
typedef struct coord            Coord;
typedef void                    Param;
typedef struct FXParam          FXParam;
typedef float                   RSParam;
typedef struct seqprof          Seqprof;
typedef struct scor_set         Scor_set;
typedef struct clssfcn          Clssfcn;
/*typedef struct RPoint     RPoint;*/

static const int PERL_OK = 1;
static const int PERL_FAIL = 0;

MODULE = Wurst        PACKAGE = Wurst
PROTOTYPES: ENABLE

double
constant(sv,arg)
    PREINIT:
    STRLEN        len;
    INPUT:
    SV *        sv
    char *        s = SvPV(sv, len);
    int        arg
    CODE:
    RETVAL = constant(s,len,arg);
    OUTPUT:
    RETVAL


#/* ---------------- File ops     ------------------------------ */
int
file_pre_read (const char *fname)
    CODE:
        RETVAL = file_pre_read (fname);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

#/* ---------------- Sec struct data --------------------------- */

Sec_s_data *
sec_s_data_read (fname)
    char *fname;

char *
sec_s_data_string (sec_s_data)
    Sec_s_data *sec_s_data;

MODULE = Wurst        PACKAGE = Sec_s_dataPtr PREFIX = sec_s_data_

void
sec_s_data_DESTROY (sec_s_data)
        Sec_s_data *sec_s_data
    CODE:
        sec_s_data_destroy (sec_s_data);

MODULE = Wurst        PACKAGE = Wurst
#/* ---------------- Sequences  -------------------------------- */
Seq *
seq_read(fname)
    char *fname;

Seq *
seq_from_string (s)
        char *s;

Seq *
seq_get_1 ( ps_a, n)
       Seq_array **ps_a
       size_t n;

size_t
seq_num ( ps_a)
        Seq_array **ps_a;

size_t
seq_size ( seq)
        Seq *seq;


char *
seq_print ( seq )
        Seq *seq

Seq_array **
seq_read_many (fname, s_a = NULL)
        const char *fname
        Seq_array **s_a;

char *
seq_print_many (ps_a)
        Seq_array **ps_a

Seq *
seq_duplicate( src, n)
        Seq *src
        size_t n

MODULE = Wurst       PACKAGE = SeqPtr  PREFIX = seq_

void
seq_DESTROY(seq)
       Seq *seq
    CODE:
        seq_destroy (seq);

MODULE = Wurst       PACKAGE = Seq_arrayPtrPtr PREFIX = seq_arrayPtr_

void
seq_arrayPtr_DESTROY(s_a)
        Seq_array **s_a
    CODE:
        seq_array_destroy (s_a);


MODULE = Wurst        PACKAGE = Wurst

#/* ---------------- sub_mat -----------------------------------
# */

Sub_mat *
sub_mat_read (fname)
    const char *fname

char *
sub_mat_string(s_m)
    Sub_mat *s_m

Sub_mat *
sub_mat_shift (s_m, x)
    Sub_mat *s_m
    float x;

void
sub_mat_scale (s_m, bot, top)
    Sub_mat *s_m
    int bot
    int top

float
sub_mat_get_by_i (s_m, ndx1, ndx2)
    Sub_mat *s_m
    int ndx1
    int ndx2

float
sub_mat_get_by_c (s_m, ndx1, ndx2)
    Sub_mat *s_m
    char ndx1
    char ndx2

int
sub_mat_set_by_i (s_m, ndx1, ndx2, f)
        Sub_mat *s_m
        int ndx1
        int ndx2
        float f
    CODE:
        RETVAL = sub_mat_set_by_i (s_m, ndx1, ndx2, f);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

int
sub_mat_set_by_c (s_m, ndx1, ndx2, f)
        Sub_mat *s_m
        char ndx1
        char ndx2
        float f
    CODE:
        RETVAL = sub_mat_set_by_c (s_m, ndx1, ndx2, f);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

MODULE = Wurst       PACKAGE = Sub_matPtr PREFIX = sub_mat_

void
sub_mat_DESTROY (s_m)
        Sub_mat *s_m
    CODE:
        sub_mat_destroy (s_m);


#/* ---------------- coord   -----------------------------------
# */

MODULE = Wurst PACKAGE = Wurst

int
coord_2_bin (c, fname)
        Coord *c;
        const char *fname;
    CODE:
        RETVAL = coord_2_bin (c, fname);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

char *
coord_name (c)
    Coord *c;

size_t
coord_size (c)
    Coord *c;

Coord *
coord_read ( fname )
    const char *fname;

Seq *
coord_get_seq (c)
        Coord *c;

int
coord_has_sec_s (c)
        Coord *c;

char *
coord_get_sec_s(c)
        Coord *c;

int
coord_2_pdb ( fname, c, seq = NULL)
        const char *fname
        Coord *c
        Seq   *seq;
    CODE:
        RETVAL = coord_2_pdb (fname, c, seq);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

int
coord_2_spdb ( fname, c, scorset, seq = NULL)
     const char *fname;
     Coord *c
     Scor_set *scorset
     Seq   *seq;
    CODE:
        RETVAL = coord_2_spdb (fname, c, seq, scorset);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

float
coord_c_n_dist (c, i, j, sqrt_flag)
        Coord *c
        unsigned int i
        unsigned int j
        unsigned int sqrt_flag


void
coord_calc_psi (c)
        Coord *c


void
coord_calc_phi (c)
        Coord *c


float
coord_psi (c, j, shift_min);
        Coord *c
        const size_t j
        const float shift_min
      CODE:
       RETVAL = coord_psi (c, j, shift_min);
       if (RETVAL == BAD_ANGLE)
           XSRETURN_UNDEF;
      OUTPUT:
        RETVAL


float
coord_phi (c, j, shift_min);
        Coord *c
        const size_t j
        const float shift_min
      CODE:
       RETVAL = coord_phi (c, j, shift_min);
       if (RETVAL == BAD_ANGLE)
           XSRETURN_UNDEF;
      OUTPUT:
        RETVAL

int
coord_geo_gap (c, scale, max)
        Coord *c
        float scale
        float max
    PPCODE:
    {
        unsigned int num_gap;
        int r;
        float quad, linear, logistic;
        SV *sv1 = &PL_sv_undef;
        SV *sv2 = &PL_sv_undef;
        SV *sv3 = &PL_sv_undef;
        SV *sv4 = &PL_sv_undef;
        r = coord_geo_gap (c, &quad, &linear, &logistic, &num_gap, scale, max);
        if (r == EXIT_FAILURE)
            XSRETURN_UNDEF;
        sv1 = newSVnv (quad);
        sv2 = newSVnv (linear);
        sv3 = newSVnv (logistic);
        sv4 = newSViv (num_gap);
        sv_2mortal (sv1);
        sv_2mortal (sv2);
        sv_2mortal (sv3);
        sv_2mortal (sv4);
        XPUSHs (sv1);
        XPUSHs (sv2);
        XPUSHs (sv3);
        XPUSHs (sv4);
    }

int
dme_thresh (frac, c1, c2, thresh)
        float frac = NO_INIT
        Coord *c1
        Coord *c2
        float thresh;
    CODE:
        RETVAL = dme_thresh (&frac, c1, c2, thresh);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;
    OUTPUT:
        frac
        RETVAL

Coord *
pdb_read ( fname, acq_c, chain)
    const char *fname
    const char *acq_c
    const char chain;


MODULE = Wurst  PACKAGE = CoordPtr  PREFIX = coord_

void
coord_DESTROY (c)
        Coord *c
    CODE:
        coord_destroy (c);

MODULE = Wurst  PACKAGE = Wurst
#/* ---------------- make_model --------------------------------
# */

Coord *
make_model (p, seq, coord)
        Pair_set *p
        Seq      *seq
        Coord    *coord


int
mc_run (fname, frag_len, max_iter)
        const char *fname
        const int frag_len
        const int max_iter

int
model_pdb_num ( mdl, resnum)
        Coord    *mdl
        int resnum;
     CODE:
        RETVAL = model_pdb_num(mdl, resnum);
        if (RETVAL == -99999)
            XSRETURN_UNDEF;
     OUTPUT:
        RETVAL

int
model_res_num ( mdl, mnum)
        Coord    *mdl
        int   mnum;
     CODE:
        RETVAL = model_res_num(mdl, mnum);
        if (RETVAL == -1)
          XSRETURN_UNDEF;

float *
coord_2_pnlty (c, f)
        Coord *c
        float f


#/* ---------------- scoring routines --------------------------
# */
MODULE = Wurst PACKAGE = Wurst

Score_mat *
score_mat_new (n_rows, n_cols)
    size_t n_rows
    size_t n_cols

Score_mat *
score_mat_read (fname)
    const char *fname

Score_mat *
score_mat_add (mat1, mat2, scale, shift = 0.0)
        Score_mat *mat1
        Score_mat *mat2
        float scale
        float shift;

void
score_mat_info (IN Score_mat *mat, \
                OUTLIST float min, OUTLIST float max, OUTLIST float av, \
                OUTLIST float std_dev)
    PROTOTYPE: $
    CODE:
        score_mat_info (mat, &min, &max, &av, &std_dev);

Score_mat *
score_mat_scale (mat, scale)
        Score_mat *mat
        const float scale;

Score_mat *
score_mat_shift (mat1, shift)
        Score_mat *mat1
        const float shift;

int
score_mat_write ( mat1, fname)
        Score_mat *mat1
        const char *fname;
    CODE:
        RETVAL = score_mat_write (mat1, fname);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

void
score_mat_diag_wipe( pair_set, smat )
            Pair_set *pair_set
            Score_mat *smat;

Score_mat *
score_mat_double_matrix ( smat )
            Score_mat *smat;


int
score_mat_write_gnuplot( smat, fname, protA, protB )
            Score_mat  *smat
            const char *fname
            const char *protA
            const char *protB;


int
score_dpt (score_mat, seq1, seq2, d_p)
        Score_mat *score_mat
        Seq *seq1
        Seq *seq2
        Dpt_list *d_p;
    CODE:
        RETVAL = score_dpt (score_mat, seq1, seq2, d_p);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

int
score_prof_prof (score_mat, sp1, sp2, subst_matrix)
        Score_mat *score_mat
        Seqprof *sp1
        Seqprof *sp2
        Sub_mat *subst_matrix;
    CODE:
        RETVAL = score_prof_prof (score_mat, sp1, sp2, subst_matrix);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

int
score_smat (score_mat, seq1, seq2, subst_matrix)
        Score_mat *score_mat
        Seq *seq1
        Seq *seq2
        Sub_mat *subst_matrix;
    CODE:
        RETVAL = score_smat (score_mat, seq1, seq2, subst_matrix);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

int
score_sprof (score_mat, sp, seq, subst_matrix)
        Score_mat *score_mat
        Seqprof *sp
        Seq *seq
        Sub_mat *subst_matrix;
    CODE:
        RETVAL = score_sprof (score_mat, sp, seq, subst_matrix);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;


int
score_fx (score_mat, seq, coord, param)
        Score_mat *score_mat
        Seq *seq
        Coord *coord
        FXParam *param
    CODE:
        RETVAL = score_fx (score_mat, seq, coord, param);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

int
score_fx_prof (score_mat, sp, coord, param)
        Score_mat *score_mat
        Seqprof *sp
        Coord *coord
        FXParam *param
    CODE:
        RETVAL = score_fx_prof (score_mat, sp, coord, param);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;


int
score_sec (score_mat, sec_s, coord)
        Score_mat *score_mat
        Sec_s_data *sec_s
        Coord *coord
    CODE:
        RETVAL = score_sec (score_mat, sec_s, coord);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

RSParam *
param_rs_read(fname)
    const char *fname;

float
score_rs (coord, params)
    Coord *coord
    RSParam *params;



FXParam *
param_fx_read (fname)
   const char *fname


MODULE = Wurst  PACKAGE = RSParamPtr  PREFIX = RSParam_
void
RSParam_DESTROY (p)
        RSParam *p
    CODE:
        param_rs_destroy (p);


MODULE = Wurst  PACKAGE = FXParamPtr  PREFIX = FXParam_

void
FXParam_DESTROY (f)
        FXParam *f
    CODE:
        FXParam_destroy(f);

MODULE = Wurst PACKAGE = Wurst

char *
score_mat_string (scores, seq0, seq1)
    Score_mat *scores
    Seq *seq0
    Seq *seq1;

Score_mat *
score_pvec_quadrifiedsmat (score_mat, pvec1, pvec2)
    Score_mat *score_mat
    Prob_vec *pvec1
    Prob_vec *pvec2;

Pair_set *
score_mat_sum_full (rmat, scores, pgap_open, pgap_widen, qgap_open, qgap_widen,p_mult, q_mult, align_type, bias_set = NULL)
        Score_mat * &rmat = NO_INIT
        Score_mat *scores
        float pgap_open
        float pgap_widen
        float qgap_open
        float qgap_widen
        float *p_mult
        float *q_mult
        int align_type
        Pair_set *bias_set
    OUTPUT:
        rmat
        RETVAL


MODULE = Wurst    PACKAGE = floatPtr     PREFIX = float_
void
float_DESTROY (f)
       float *f
    CODE:
       if (f)
           free (f);

MODULE = Wurst    PACKAGE = Score_matPtr PREFIX = score_mat_
void
score_mat_DESTROY (s)
        Score_mat *s
    CODE:
        score_mat_destroy (s);

#/* ---------------- pair_set ---------------------------------
#/* This is very messy. We want to return a list with two strings
#/* on it. Current versions of xsubpp let us use the nicer IN/OUT
#/* keywords. Older versions require us to PPCODE and dork around
#/* in the return stack explicitly.
# */
MODULE = Wurst PACKAGE = Wurst

char *
pair_set_chimera (s, c1, c2)
        Pair_set *s
        Coord *c1
        Coord *c2


void
pair_set_coverage (s, n1, n2)
        Pair_set *s
        size_t n1
        size_t n2
    PPCODE:
    {
        char *c1, *c2;
        if (pair_set_coverage (s, n1, n2, &c1, &c2) ==  EXIT_FAILURE) {
            free (c1);
            free (c2);
            XSRETURN_UNDEF;
        }
        EXTEND(SP, 2);
        PUSHs(sv_2mortal(newSVpv (c1, 0)));
        PUSHs(sv_2mortal(newSVpv (c2, 0)));
        free (c1);
        free (c2);
    }

int
pair_set_gap (s, open_scale, widen_scale)
        Pair_set *s
        float open_scale
        float widen_scale
    PPCODE:
    {
        float open_cost, widen_cost;
        int r;
        r = pair_set_gap (s, &open_cost, &widen_cost, open_scale, widen_scale);
        if (r == EXIT_FAILURE)
            XSRETURN_UNDEF;
        XPUSHs(sv_2mortal(newSVnv (open_cost)));
        XPUSHs(sv_2mortal(newSVnv (widen_cost)));
    }

char *
pair_set_string (s, seq0, seq1)
        Pair_set *s
        Seq *seq0
        Seq *seq1;

char *
multal_string (s)
        Pair_set *s;

char *
pair_set_pretty_string (s, seq0, seq1, sec_s_data = NULL, c2 = NULL)
        Pair_set *s
        Seq *seq0
        Seq *seq1
        Sec_s_data *sec_s_data
        Coord *c2;

void
pair_set_get_alignment_indices (IN Pair_set *p_s, IN int sequencenumber,\
                              OUTLIST int start, OUTLIST int stop)
         PROTOTYPE: $$
         CODE:
                   pair_set_get_alignment_indices (p_s, sequencenumber, &start, &stop);

int
pair_set_circularpermutated( p_s , prot1_len , prot2_len )
         Pair_set *p_s
         size_t   prot1_len
         size_t   prot2_len;


int
get_seq_id_simple (  Pair_set *pair_set, Seq *s1, Seq *s2 );

#/* As above, we would like to use the newer IN/OUT lists, but
#/* for compatibility with older perl releases, we use PPCODE and take
#/* care ofthe stack explicitly
# */
int
pair_set_score(s)
        Pair_set *s;
    PPCODE:
    {
        float score, score_smpl;
        SV *sv1 = &PL_sv_undef;
        SV *sv2 = &PL_sv_undef;
        if (pair_set_score (s, &score, &score_smpl) == EXIT_FAILURE)
            XSRETURN_UNDEF;
        sv1 = newSVnv (score);
        sv_2mortal (sv1);
        sv2 = newSVnv (score_smpl);
        sv_2mortal (sv2);
        XPUSHs (sv1);
        XPUSHs (sv2);
        RETVAL = 1;
    }

int
pair_set_extend(s, n0, n1, ext_len)
        Pair_set *s
        size_t n0
        size_t n1
        int ext_len;

MODULE = Wurst  PACKAGE = Pair_setPtr  PREFIX = pair_set_

void
pair_set_DESTROY (s)
        Pair_set *s
    CODE:
        pair_set_destroy (s);

MODULE = Wurst PACKAGE = Wurst
int
func_int ( )

float
func_float ()

char *
func_char ()

char *
funcs1_char (in)
        char *in;

char *
funcs2_char ()

#/* ---------------- clean_up   --------------------------------
void
free_scratch ()


#/* ---------------- temporary solvation    --------------------
# * We have some very temporary solvation routines. When we have
# * decided on an interface, we should put them in properly.
# * For the moment, these things just write to their own files
# * in a rather uncontrolled manner.
# */

int
get_nbor (c)
    Coord * c

int
get_solv_data (c)
	Coord *c


int
get_torsion_data (c)
	Coord *c


#/* ---------------- blast checkpoint files --------------------
# */

Seqprof *
blst_chk_read (f)
        const char *f;

char *
seqprof_str (p)
        Seqprof *p;

Seq *
seqprof_get_seq (c)
        Seqprof *c;


MODULE = Wurst PACKAGE = SeqprofPtr PREFIX = seqprof_
void
seqprof_DESTROY (p)
        Seqprof *p
    CODE:
        seqprof_destroy (p);

#/* -------------------- scor_set routines ---------------
# */

MODULE = Wurst PACKAGE = Wurst

void
scor_set_fromvec ( dblvec )
    AV *dblvec;
    PREINIT:
       size_t s=0;
       size_t t;
       SV **d_el, *ssv;
       Scor_set *ss=NULL;
       double *dv;
    PPCODE:
       t = 1+av_len(dblvec);
       if (t>0) {
         dv = malloc(sizeof(*dv)*t); /* vector passed to C function */
         while (s<t) {
           d_el=av_fetch(dblvec, s, 0);
           if (d_el && (SvNIOK(*d_el)))
             /* copy the values from the perl array */
             dv[s] = SvNV( *d_el );
           else
             dv[s] = 0.0; /* and do not complain if */
                         /* we're given garbage   */
           s++;
         }
         ss = scor_set_fromvec( t, dv);
         free(dv);
       }
       if (ss) {
           ssv = sv_newmortal();
           sv_setref_pv(ssv, "Scor_setPtr", (void *)ss);
           EXTEND(SP,1);
           PUSHs(ssv);
       } else {
           XSRETURN_UNDEF;
       }


MODULE = Wurst PACKAGE = Wurst

AV*
scor_set_to_vec ( scoreset )
       Scor_set *scoreset;

       PREINIT:
          AV *av=Nullav;
          SV *a_el;
          size_t i=0;

       CODE:
          if (scoreset!=NULL) {
             av = newAV();
             av_unshift(av, scoreset->n);
             while (i<scoreset->n) {
                a_el = sv_newref(newSVnv(scoreset->scores[i]));
                if (av_store(av, i, a_el)==NULL)
                     croak("scor_set_to_vec: av_store returned null - \
 where is my array ???");

                i++;
             }
           }
           RETVAL = av;
        OUTPUT:
           RETVAL


Scor_set *
scor_set_simpl ( pairset, scoremat)
        Pair_set *pairset
        Score_mat *scoremat;

int
scor_set_scale  ( scorset, scale )
        Scor_set *scorset
        float     scale;
    CODE:
        RETVAL = scor_set_scale( scorset, scale);
        if (RETVAL == EXIT_FAILURE)
          XSRETURN_UNDEF;
        else
          RETVAL = 1;


MODULE = Wurst PACKAGE = Scor_setPtr PREFIX = scor_set_

void
scor_set_DESTROY (p)
        Scor_set *p
    CODE:
        scor_set_destroy (p);

#/* ---------------- class probability vectors -----------------
# */
MODULE = Wurst PACKAGE = Wurst

char *
prob_vec_info (pvec)
    Prob_vec *pvec;

#Prob_vec *
#seq_2_prob_vec ( seq, ac )
#    Seq *seq
#    Aa_clssfcn *ac;

int
score_pvec (score_mat, pvec1, pvec2)
    Score_mat *score_mat
    Prob_vec *pvec1
    Prob_vec *pvec2
    CODE:
        RETVAL = score_pvec (score_mat, pvec1, pvec2);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

size_t
prob_vec_size(pvec);
    Prob_vec * pvec;

size_t
prob_vec_length(pvec);
    Prob_vec * pvec;

Prob_vec *
prob_vec_copy(pvec);
    Prob_vec * pvec;

int
prob_vec_write (pvec, fname)
    Prob_vec *pvec
    const char *fname
    CODE:
        RETVAL = prob_vec_write (pvec, fname);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

Prob_vec *
prob_vec_read (fname)
    const char *fname

Prob_vec *
prob_vec_duplicate( pvec, k);
        Prob_vec * pvec
        size_t k;

MODULE = Wurst   PACKAGE = Prob_vecPtr PREFIX = prob_vec_
void
prob_vec_DESTROY (pvec)
        Prob_vec *pvec
    CODE:
        prob_vec_destroy (pvec);

#/* ---------------- amino acid + structure classification -----
# */
MODULE=Wurst PACKAGE = Wurst

Aa_strct_clssfcn *
aa_strct_clssfcn_read (fname, error)
    const char *fname
    const float error

void
aa_strct_dump (ac)
    Aa_strct_clssfcn *ac

size_t
aa_strct_nclass (ac)
    Aa_strct_clssfcn *ac

size_t
aa_strct_size (ac)
    Aa_strct_clssfcn *ac

Prob_vec *
strct_2_prob_vec (structure, cmodel, norm)
    Coord * structure
    Aa_strct_clssfcn * cmodel
    int norm

Prob_vec *
strct_2_duplicated_prob_vec ( structure, seq, sp, size, cmodel, n_duplications)
    Coord * structure
    Seq * seq
    Seqprof * sp
    size_t size
    Aa_strct_clssfcn * cmodel
    size_t n_duplications;

Prob_vec *
aa_strct_2_prob_vec (structure, cmodel, norm)
    Coord * structure
    Aa_strct_clssfcn * cmodel
    int norm

Prob_vec *
prof_aa_strct_2_prob_vec (structure, sp, cmodel, norm)
    Coord * structure
    Seqprof * sp
    Aa_strct_clssfcn * cmodel
    int norm

Prob_vec *
aa_2_prob_vec (seq, cmodel, norm)
    Seq * seq
    Aa_strct_clssfcn * cmodel
    int norm

Prob_vec *
prof_aa_2_prob_vec (sp, cmodel, norm);
    Seqprof * sp
    Aa_strct_clssfcn * cmodel
    int norm;


#MODULE = Wurst   PACKAGE = Aa_strct_clssfcnPtr PREFIX = aa_strct_clssfcn_

void
aa_strct_clssfcn_DESTROY (aa)
        Aa_strct_clssfcn *aa
    CODE:
        aa_strct_clssfcn_destroy (aa);


#/* ---------------- amino acid classification------------------
# */
MODULE = Wurst PACKAGE = Wurst

Aa_clssfcn *
ac_read (fname)
    const char *fname

size_t
ac_size (ac)
    Aa_clssfcn *ac

size_t
ac_nclass (ac)
    Aa_clssfcn *ac

void
ac_dump (ac)
    Aa_clssfcn *ac

int
computeMembershipAAProf (mship, sp, ac)
        float ** mship
        Seqprof * sp
        Aa_clssfcn * ac
    CODE:
        RETVAL = computeMembershipAAProf (mship, sp, ac);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

int
computeMembershipAA (mship, seq, ac)
    float ** mship
    Seq * seq
    Aa_clssfcn * ac
    CODE:
        RETVAL = computeMembershipAA (mship, seq, ac);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

Prob_vec *
seq_2_prob_vec (seq, ac)
    Seq *seq
    Aa_clssfcn *ac

MODULE = Wurst   PACKAGE = Aa_clssfcnPtr PREFIX = aa_clssfcn_
void
aa_clssfcn_DESTROY (aa)
        Aa_clssfcn *aa
    CODE:
        aa_clssfcn_destroy (aa);

#/* ---------------- dipeptide routines     --------------------
# */


MODULE = Wurst PACKAGE = Wurst
Dpt_list *
dpt_read   ( fname )
    const char *fname


char *
dpt_string (d_p)
         Dpt_list *d_p;

int
dpt_get_n (d_p)
         Dpt_list *d_p;

float
dpt_get_val (d_p, n)
         Dpt_list *d_p
         size_t n;
    CODE:
        {
            int error;
            RETVAL = dpt_get_val (d_p, n, &error);
            if (error == EXIT_FAILURE)
                XSRETURN_UNDEF;
        }
    OUTPUT:
        RETVAL

int
dpt_set_val (d_p, n, val)
        Dpt_list *d_p
        size_t n
        float val
    CODE:
        RETVAL = dpt_set_val (d_p, n, val);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
        else
            RETVAL = 1;

MODULE = Wurst   PACKAGE = Dpt_listPtr PREFIX = dpt_list_
void
dpt_list_DESTROY (d_p)
        Dpt_list *d_p
    CODE:
        dpt_list_destroy (d_p);


#/* ---------------- alternative scores routines  --------------------
# */
MODULE = Wurst PACKAGE = Wurst


float
find_alt_path_score_simple (Score_mat *score_mat, Pair_set *pair_set);


#/* ---------------- bayesian fragments routines  --------------------
# */
MODULE = Wurst PACKAGE = Wurst

#Prob_vec *
#struct_2_prob_vec(structure, cmodel);
#    Coord * structure
#    Clssfcn * cmodel;

Clssfcn *
get_clssfcn(influence_report_filename, abs_error);
    const char * influence_report_filename
    float abs_error;

float *
computeMembership(mship, test_vec, cmodel);
    float * mship
    float * test_vec
    Clssfcn * cmodel;


MODULE = Wurst   PACKAGE = ClssfcnPtr PREFIX = clssfcn_
void
clssfcn_DESTROY (c)
        Clssfcn *c
    CODE:
        clssfcn_destroy (c);

#/* ---------------- multiple alignments ------------------
# */
MODULE = Wurst PACKAGE = Wurst

Pair_set *
merge_alignments (align1, align2, alignment);
        Pair_set * align1
        Pair_set * align2
        Pair_set * alignment;


Prob_vec *
pvec_avg (p_vec1, p_vec2, p_set, cur_step);
        Prob_vec * p_vec1
        Prob_vec * p_vec2
		Pair_set * p_set
		int cur_step;

Pair_set *
remove_seq (pairset, index);
		Pair_set * pairset
		int index;

Pair_set *
split_multal (pairset, a, b);
		Pair_set * pairset
		int a
		int b;

#/* ---------------- superimposing structures, RMSD ------------------
# */
MODULE = Wurst PACKAGE = Wurst

NO_OUTPUT int
coord_rmsd ( IN Pair_set *pair_set, IN Coord *coord1, IN Coord *coord2, \
            IN int sub_flag,\
            OUTLIST float rmsd, OUTLIST Coord *c1_new, OUTLIST Coord *c2_new)
    PROTOTYPE: $$$$
    CODE:
        RETVAL = coord_rmsd(pair_set, coord1, coord2, sub_flag, &rmsd, &c1_new, &c2_new);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;

NO_OUTPUT int 
get_rmsd(IN Pair_set *pairset, IN Coord *r1, \
             IN Coord *r2, OUTLIST float rmsd, OUTLIST int count)
	PROTOTYPE: $$$
	CODE:
		RETVAL = get_rmsd(pairset, r1, r2, &rmsd, &count);
        if (RETVAL == EXIT_FAILURE)
            XSRETURN_UNDEF;
