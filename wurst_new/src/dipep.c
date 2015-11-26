/*
 * 27 Jan 2005
 * Read the information in a dipeptide comparison file.
 * This is a substitution matrix for dipeptides. Right now, it is
 * sparse, so we only read up the few elements that are set.
 * Internally, we store data using our own labels ("thomas"
 * format) for amino acids.
 */

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "amino_a.h"
#include "dipep_i.h"
#include "e_malloc.h"
#include "fio.h"
#include "misc.h"
#include "mprintf.h"
#include "read_seq_i.h"
#include "scratch.h"


#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: dipep.c,v 1.1 2005/03/02 10:06:39 torda Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */
 
/* ---------------- structures -------------------------------- */
struct dipep {
    float val;        /* similarity*/
    unsigned short k; /* separation of amino acids. k=1 is adjacent */
    char a1;
    char a2;
    char b1;
    char b2;
};
struct dpt_list {
    size_t n;
    struct dipep *dipep;
};

struct two_aa {
    char a;
    char b;
};

/* ---------------- valid_dipep -------------------------------
 * Check if a string is of the form "AB" or "ab" or something
 * that looks like a reasonable dipeptide.
 * Next, check if the amino acids are valid.
 * Convert the names to our internal (so-called "thomas") format
 * and put the converted amino acids into the two_aa structure.
 */
static int
valid_dipep (const char *s, struct two_aa *two_aa)
{
    const char *this_sub = "valid_dipep";
    const char *invalid_aa = "%c is not a valid amino acid\n";
    if (strlen (s) != 2) {
        err_printf (this_sub, "String %s should be 2 char long\n", s);
        return EXIT_FAILURE;
    }
    if (aa_invalid (*s)) {
        err_printf (this_sub, invalid_aa, *s);
        return EXIT_FAILURE;
    }
    if (aa_invalid (*(s+1))) {
        err_printf (this_sub, invalid_aa, *s);
        return EXIT_FAILURE;
    }
    two_aa->a = std2thomas_char (*s);
    two_aa->b = std2thomas_char (*(s+1));
    return EXIT_SUCCESS;
}

/* ---------------- add_dpt   ---------------------------------
 * The format is a whitespace-separated line like
 * ac de 5 4.5
 * Most of the code is spent checking if characters are valid. If
 * so, we build up a single piece of data in temp_dipep.
 * If we get to the end without any errors, we take the existing
 * list, call realloc to extend it and copy over the new entry.
 */
static int
add_dpt (struct dpt_list *dpt_list, char inbuf[],
         int nr_line, const char *fname)
{
    struct two_aa two_aa;
    long int tmp_l;
    size_t to_mall;
    const char *this_sub = "add_dpt";
    const char *sep = " 	";
    char *s, *t, *dist, *val;
    
    
    struct dipep temp_dipep;
    memset (&temp_dipep, 0, sizeof (temp_dipep));
    if ((s = strtok (inbuf, sep)) == NULL)
        goto error;
    if ((t = strtok (NULL, sep)) == NULL)
        goto error;
    if ((dist = strtok (NULL, sep)) == NULL)
        goto error;
    val = strtok (NULL, sep);           /* no error if last value is missing */
    if (valid_dipep ( s, &two_aa) == EXIT_FAILURE)
        goto error;
    temp_dipep.a1 = two_aa.a;
    temp_dipep.a2 = two_aa.b;
    if (valid_dipep ( t, &two_aa) == EXIT_FAILURE)
        goto error;
    temp_dipep.b1 = two_aa.a;
    temp_dipep.b2 = two_aa.b;
    tmp_l = strtol(dist, NULL, 10);
    if (tmp_l > USHRT_MAX || tmp_l < 1) {
        err_printf (this_sub, "separation \"%s\" should be small +ve\n", dist);
        goto error;
    }
    temp_dipep.k = (unsigned short) tmp_l;
    if (val) {
        double tmp;
        errno = 0;
        tmp = strtod (val, NULL);
        if (errno && (tmp == 0)) {
            err_printf (this_sub, "invalid float number %s\n", val);
            goto error;
        }
        temp_dipep.val = (float) tmp;
    } else {
        temp_dipep.val = 1.0;
    }

    dpt_list->n++;
    to_mall = (dpt_list->n * sizeof (temp_dipep));
    dpt_list->dipep = E_REALLOC (dpt_list->dipep, to_mall);
    dpt_list->dipep[ dpt_list->n-1 ] = temp_dipep;
    return EXIT_SUCCESS;
 error:
    err_printf (this_sub, "Bad input file \"%s\": %d\n", fname, nr_line);
    return EXIT_FAILURE;

}

/* ---------------- dpt_cmpr   --------------------------------
 * We have to sort our list of dipeptide information according to
 * the separation of residues in the dipeptide.
 * It seems to be legal to declare this static.
 */
static int
dpt_cmpr (const void *s1, const void *s2)
{
    const struct dipep *p1 = s1;
    const struct dipep *p2 = s2;
    const unsigned short k1 = p1->k;
    const unsigned short k2 = p2->k;
    if (k1 == k2)
        return 0;
    if (k1 < k2)
        return -1;
    return 1;
}


/* ---------------- dpt_read   --------------------------------
 * Given a filename containing dipeptide information, return
 * either a structure containing the information or a NULL
 * pointer if an error occurs.
 */
struct dpt_list *
dpt_read (const char *fname)
{
    FILE *fp;
    struct dpt_list *dpt_list = NULL;
#   ifndef BUFSIZ
        enum { BUFSIZ = 1024 };
#   endif
    char inbuf [BUFSIZ];
    int nr_line = 0;
    const char *this_sub = "dpt_read";

    if ((fp = mfopen (fname, "r", this_sub)) == NULL)
        return NULL;
    dpt_list = E_MALLOC (sizeof (*dpt_list));
    memset (dpt_list, 0, sizeof (*dpt_list));
    while (get_nline (fp, inbuf, &nr_line, BUFSIZ)) {
        if (add_dpt (dpt_list, inbuf, nr_line, fname) == EXIT_FAILURE) {
            dpt_list_destroy (dpt_list);
            return NULL;
        }
    }

    qsort (dpt_list->dipep, dpt_list->n, sizeof dpt_list->dipep[0], dpt_cmpr);
    if (dpt_list->dipep->k == 0) {
        err_printf (this_sub, "bad entry on first line of %s\n", fname);
        dpt_list_destroy (dpt_list);
        return NULL;
    }
    
    if (dpt_list->n == 0) {
        err_printf (this_sub, "No data in \"%s\"\n", fname);
        dpt_list_destroy (dpt_list);
        return NULL;
    }
    return dpt_list;
}

/* ---------------- dpt_get_n   -------------------------------
 * Given a dpt_list, return the number of elements.
 */
int
dpt_get_n (const struct dpt_list *dpt_list)
{
    return dpt_list->n;
}

/* ---------------- dpt_get_val -------------------------------
 * Return the value associated with the Nth entry on the dpt
 * list. We have to be able to return errors, so we do this via
 * the int pointer.
 */
float
dpt_get_val ( const struct dpt_list *dpt_list, size_t n, int *error)
{
    if (n >= dpt_list->n) {
        *error = EXIT_FAILURE;
        return 0;
    }
    *error = EXIT_SUCCESS;
    return dpt_list->dipep[n].val;
}

/* ---------------- dpt_set_val -------------------------------
 * Set the value associated with the Nth entry on the
 * dpt_list. Return EXIT_SUCCESS/FAILURE.
 */
int
dpt_set_val ( const struct dpt_list *dpt_list, const size_t n, const float val)
{
    if (n >= dpt_list->n)
        return EXIT_FAILURE;
    dpt_list->dipep[n].val = val;
    return EXIT_SUCCESS;
}

/* ---------------- dpt_string  -------------------------------
 * Given a dpt_list, return its contents. We use the mechanism
 * have for the string scratch buffer.
 */
char *
dpt_string ( const struct dpt_list *dpt_list)
{
    char *s;
    struct dipep *i, *ilast;
    const char *fmt = "%5c %c %5c %c  %4hu %7.5g\n";
    scr_reset();
    s = scr_printf ("%6s %7s %5s %7s\n", "#dipep1", "dipep2", "dist", "value");
    i = dpt_list->dipep;
    ilast = i + dpt_list->n;
    for ( ; i < ilast; i++) {
        char a1 = thomas2std_char (i->a1);
        char a2 = thomas2std_char (i->a2);
        char b1 = thomas2std_char (i->b1);
        char b2 = thomas2std_char (i->b2);
        s = scr_printf (fmt, a1, a2, b1, b2, i->k, i->val);
    }
    return s;
}

/* ---------------- dpt_destroy -------------------------------
 * Clean up anything associated with a dipeptide list.
 * Although this can be called from the C level, it is most
 * likely called from the interpreter when it decides to get
 * rid of a dpt_list structure.
 */
void
dpt_list_destroy (struct dpt_list *dpt_list)
{
    if (dpt_list == NULL)
        return;
    if (dpt_list->dipep != NULL)
        free (dpt_list->dipep);
    free (dpt_list);
}

/* ---------------- Start of scoring section ------------------
 * Maybe this should go into a separate file. If so, we do so
 * later. Above is the code for manipulating the dipeptide data.
 * Here is the code for doing a score matrix. In anticipation, we
 * will even have a separate #include section below.
 */

#include "score_mat.h"
#include "seq.h"


/* ---------------- print_dpt   -------------------------------
 * This is for debugging and can be removed a bit later.
 * Given a single piece of dipeptide information, print it to
 * stdout.
 */
/* #define want_print_dpt */
#ifdef  want_print_dpt
static void
print_dpt ( const struct dipep *dipep)
{
    char a1 = thomas2std_char (dipep->a1);
    char a2 = thomas2std_char (dipep->a2);
    char b1 = thomas2std_char (dipep->b1);
    char b2 = thomas2std_char (dipep->b2);
    const char *fmt = "%4c %c %5c %c  %4hu %7.5g at %4hu\n";
    mprintf (fmt, a1, a2, b1, b2, dipep->k, dipep->val, dipep->k);
}
#endif /* want_print_dpt */

/* ---------------- find_one    -------------------------------
 * For a single dipeptide, mark any relevant residues
 */
static unsigned int
find_one (const struct seq *seq, unsigned int *mark,
            const char a1, const char a2,
            const short unsigned int k)
{
    char *s  = seq->seq;
    const size_t len = seq->length - k;
    unsigned int i, count;
    count = 0;
    for (i = 0; i < len; i++)
        if (s[i] == a1)                         /* first char of dipeptide */
            if (s [i + k] == a2)                /* second char of dipeptide */
                mark[count++] = i;
    return count;
}
void
dump_f_matrix (const float **mat, const size_t n_rows, const size_t n_cols);

/* ---------------- inner_score -------------------------------
 * We do the scoring for each dipeptide separation, k,
 * separately. 
 * Most sanity checks have been done by the caller.
 * Remember, the first and last row and column are left at zero,
 * so we always have to add 1 when putting number in the matrix.
 */
static void
inner_score (float **mat, const struct seq *seq1, const struct seq *seq2,
             const struct dipep *d_ini, const struct dipep *d_end,
             const short unsigned int k)
{
    unsigned int *mark_s1_d1 = NULL,
                 *mark_s1_d2 = NULL,
                 *mark_s2_d1 = NULL,
                 *mark_s2_d2 = NULL;
    const struct dipep *d;

    if ((seq1->length <= (size_t)(k+1) ) || seq2->length <= (size_t)(k+1))
        return;
    mark_s1_d1 = E_MALLOC ((seq1->length - k) * sizeof (mark_s1_d1));
    mark_s1_d2 = E_MALLOC ((seq1->length - k) * sizeof (mark_s1_d1));
    mark_s2_d1 = E_MALLOC ((seq2->length - k) * sizeof (mark_s1_d1));
    mark_s2_d2 = E_MALLOC ((seq2->length - k) * sizeof (mark_s1_d1));
    for (d = d_ini; d < d_end; d++) {            /* loop over dipeptide info */
        unsigned int n_1_1, n_1_2, n_2_1, n_2_2;
        n_1_1 = find_one (seq1, mark_s1_d1, d->a1, d->a2, k);
        n_1_2 = find_one (seq1, mark_s1_d2, d->b1, d->b2, k);
        n_2_1 = find_one (seq2, mark_s2_d1, d->a1, d->a2, k);
        n_2_2 = find_one (seq2, mark_s2_d2, d->b1, d->b2, k);
        if (n_1_1 && n_2_2) {  /* seq 1, dipep first with seq 2 dipep second */
            unsigned short i, j;
            for (i = 0; i < n_1_1; i++) {
                unsigned short i_p = mark_s1_d1[i] + 1;   /* The +1 accounts */
                for (j = 0; j < n_2_2; j++) {             /* for the zero at */
                    unsigned short j_p = mark_s2_d2[j] + 1;  /* start of mat */
                    mat [i_p] [j_p]        = d->val;
                    mat [i_p + k][j_p + k] = d->val;
                }
            }
        }
        if (n_1_2 && n_2_1) {
            unsigned short i, j;
            for (i = 0; i < n_1_2; i++) {
                unsigned short i_p = mark_s1_d2[i] + 1;
                for (j = 0; j < n_2_1; j++) {          
                    unsigned short j_p = mark_s2_d1[j] + 1;
                    mat [i_p] [j_p]        = d->val;
                    mat [i_p + k][j_p + k] = d->val;
                }
            }
        }
    }
    free (mark_s1_d1);
    free (mark_s1_d2);
    free (mark_s2_d1);
    free (mark_s2_d2);
}

/* ---------------- score_dpt   -------------------------------
 * The first version of a score function based on dipeptides.
 * The interface should be as similar as possible to the other
 * score functions. This means we return EXIT_SUCCESS/FAILURE
 * We can assume the list of information is already sorted.
 */
int
score_dpt ( struct score_mat *score_mat, struct seq *s1,
            struct seq *s2, const struct dpt_list *dpt_list)
{
    const char *this_sub = "score_dpt";
    extern const char *mismatch;

    if (( score_mat->n_rows != s1->length + 2) ||
        ( score_mat->n_cols != s2->length + 2)) {
        err_printf (this_sub, mismatch,
                    score_mat->n_rows - 2, score_mat->n_cols - 2,
                    s1->length, s2->length);
        return EXIT_FAILURE;
    }
    if (dpt_list->n == 0) {
        err_printf (this_sub, "Warning: empty dipeptide list\n");
        return EXIT_SUCCESS;
    }

    {
        const char *too_long = "seq length %ud but max is %ud\n";
        extern const char *prog_bug;
        if (s1->length >= USHRT_MAX) {
            err_printf (this_sub, prog_bug, __FILE__, __LINE__);
            err_printf (this_sub, too_long, s1->length, USHRT_MAX);
            return (EXIT_FAILURE);
        }
        if (s2->length >= USHRT_MAX) {
            err_printf (this_sub, prog_bug, __FILE__, __LINE__);
            err_printf (this_sub, too_long, s2->length, USHRT_MAX);
            return (EXIT_FAILURE);
        }
    }

    seq_std2thomas (s1);                       /* Force both sequences into */
    seq_std2thomas (s2);                       /* our "thomas" label/format. */

    {
        struct dipep *d, *dlast, *d_old;
        short unsigned k_old = dpt_list->dipep->k;
        d_old = d = dpt_list->dipep;
        dlast = d + dpt_list->n;
        for ( ; d < dlast; d++){
            if (d->k != k_old) {
                inner_score (score_mat->mat, s1, s2, d_old, d, k_old);
                k_old = d->k;
                d_old = d;
            }
        }
        inner_score (score_mat->mat, s1, s2, d_old, dlast, k_old);
    }
    return (EXIT_SUCCESS);
}
