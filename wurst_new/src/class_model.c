#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "gsldir/gsl_sf_erf.h"
#include "gsldir/gsl_linalg.h"
#include "gsldir/gsl_blas.h"

#include "amino_a.h"
#include "bad_angle.h"
#include "e_malloc.h"
#include "mprintf.h"
#include "fio.h"
#include "bad_angle.h"
#include "class_model.h"

/* ---------------- compare -----------------------------------
 * compares to size_t values ( used for qsort() )
 */
static int
compare (const void * a, const void * b)
{
    return ( (int)(*(size_t *)a - *(size_t *)b ));
}


/* ---------------- clssfcn_destroy ----------------------------
 * frees up the memory used by some classification
 */
void
clssfcn_destroy(struct clssfcn * c)
{
    size_t i, j;
    if (c != NULL) {

        if (c->classmodel != NULL) {
            for (i = 0; i < c->n_class; i++)
                free_if_not_null(c->classmodel[i]);
            free(c->classmodel);
        }

        if (c->cov_matrix != NULL) {
            for (i = 0; i < c->n_class; i++)
                free_if_not_null(c->cov_matrix[i]);
            free(c->cov_matrix);
        }

        if (c->param != NULL) {
            for (i = 0; i < c->n_class; i++) {
                if (c->param[i] != NULL) {
                    for (j = 0; j < c->dim; j++)
                        free_if_not_null(c->param[i][j]);
                    free(c->param[i]);
                }
            }
            free(c->param);
        }

        free_if_not_null(c->class_weight);

        free(c);
    }
}

static void
scanf_error (struct clssfcn *cmodel, size_t i, size_t *attr_idx,
             const char *influence_report_filename)
{
    const char this_sub[] = "scanf_error";

    long i_ = 0;

    err_printf(this_sub,"unexpected EOF reached in file %s!\n",
               influence_report_filename);
    for (i_ = i; i_ >= 0; i_--) {
        free_if_not_null(cmodel->classmodel[i_]);
        free_if_not_null(cmodel->param[i_]);
    }
    free_if_not_null(cmodel->classmodel);
    free_if_not_null(cmodel->param);
    free_if_not_null(attr_idx);
    free_if_not_null(cmodel->class_weight);
    free_if_not_null(cmodel);
}

/* ---------------- get_classfcn -------------------------------
 * Reads an influence report (only report_mode = "data" generated by
 * AutoClass-C
 * Do not forget to give an absolute error value from your original data
 * It is used for the integration in computeMembership().
 */
struct clssfcn *
get_clssfcn(const char *influence_report_filename, const float abs_error)
{
    const char this_sub[] = "get_clssfcn";
    size_t dummy_i = 0;
    float dummy_f = 0.0;
    char dummy_c;
    char dummy_s[5];
    char dummy_str[] = {'\0','\0','\0','\0','\0','\0','\0','\0','\0','\0'};
    float mean = 0, dev = 0;
    size_t corr_idx = 0;
    size_t i = 0, i_=0, j = 0, j_ = 0, j__ = 0, k = 0, l = 0;
    fpos_t pos, pos_, pos__;
    
    
    FILE * influence_report_file = NULL;
    
    struct clssfcn *cmodel = E_CALLOC( 1, sizeof(struct clssfcn));
    
    size_t n_cases = 0; /* (used to calculate a normalized class weight) */
    size_t n_attr = 0;       /* mapping between attribute indices and the  */
    size_t *attr_idx = NULL; /* actual array indices for storing in cmodel */
    int *use_idx = NULL; /* which part of the correlation matrix is found? */

    if ((influence_report_file = mfopen(influence_report_filename, "r",
                                        this_sub)) == NULL) {
        return NULL;
    }

    cmodel->abs_error = (float) fabs(abs_error);
    
    /* looks for the number of cases */
    while (fscanf(influence_report_file," DATA_CLSF_HEADER %*s %*s %*s %*s %u",
                  (unsigned int*)&n_cases) != 1) {
        if (feof(influence_report_file)) {
            err_printf(this_sub,"unexpected EOF reached in file %s!\n",
                       influence_report_filename);
            free_if_not_null(cmodel);
            return NULL;
        }
    }

    /* looks for the number of classes */
    while (fscanf(influence_report_file,
                  "%*s DATA_POP_CLASSES %*s %*s %u",
                  (unsigned int*)&cmodel->n_class) != 1) {
        if (feof(influence_report_file)) {
            err_printf(this_sub,"unexpected EOF reached in file %s!\n",
                       influence_report_filename);
            free_if_not_null(cmodel);
            return NULL;
        }
    }

    if ( cmodel->n_class < 1 ) {
        err_printf(this_sub,"not even a single class found in file %s!\n",
                   influence_report_filename);
        free_if_not_null(cmodel);
        return NULL;
    }


    cmodel->class_weight = E_CALLOC( cmodel->n_class, sizeof(float));
    cmodel->cov_matrix = E_CALLOC( cmodel->n_class, sizeof(float*));
    
    /* looks for the beginning of attribute list */
    while (fscanf(influence_report_file,
                  "%*s DATA_NORM_INF_VALS "
                  "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s "
                  "num description %s",
                  dummy_str) != 1 && dummy_str != "I-*k") {
        if (feof(influence_report_file)) {
            err_printf(this_sub,"unexpected EOF reached in file %s!\n",
                       influence_report_filename);
            free_if_not_null(cmodel->class_weight);
            free_if_not_null(cmodel);
            return NULL;
        }
        
    }

    /* Counts the number of actually used attributes, that is the dimension */
    fgetpos(influence_report_file, &pos);
    while (fscanf(influence_report_file, "%*u %5c %f", dummy_s, &dummy_f) == 2)
        if (dummy_s[0] == 'p')
            cmodel->dim++;
    
    if (fsetpos(influence_report_file, &pos) == -1) {
        err_printf (this_sub, "Fail in fsetpos in influence file\n");
        fclose (influence_report_file);
        free (cmodel);
        return NULL;
    }

    n_attr = cmodel->dim;

    use_idx = E_CALLOC( cmodel->dim, sizeof(int));
    attr_idx = E_CALLOC( n_attr, sizeof(size_t));


    /* Reads the list of attributes and sorts it to get a mapping between
     * attribute indices and the array to store them */
    i = 0;

    while (fscanf(influence_report_file, "%u %5c %f",
                  (unsigned int*)&j, dummy_s, &dummy_f) == 3) {
        if (dummy_s[0] == 'p'){
            attr_idx[i++] = j;
            fgetpos(influence_report_file, &pos_);
        }
    }

    fsetpos(influence_report_file, &pos_);

    qsort(attr_idx, n_attr, sizeof(attr_idx[0]), compare);

    cmodel->param = E_CALLOC( cmodel->n_class, sizeof(double**));

    cmodel->classmodel = E_CALLOC( cmodel->n_class, sizeof(int*));

    /* Reads the descriptions for all classes */
    for (i = 0; i < cmodel->n_class; i++) {
        cmodel->param[i] = E_CALLOC( cmodel->dim, sizeof(double*));
        cmodel->cov_matrix[i] = E_CALLOC( cmodel->dim*cmodel->dim,
                                          sizeof(double));

        for (j = 0; j < cmodel->dim; j++)
            cmodel->param[i][j] = E_CALLOC( 3, sizeof(double));

        cmodel->classmodel[i] = E_CALLOC( cmodel->dim, sizeof(int));


        /* finds the beginning of a class discription and stores its weight */
        while (fscanf(influence_report_file,
                      "%s DATA_CLASS %*u CLASS %*u - weight %f",
                      dummy_str, &cmodel->class_weight[i]) != 2) {
            if (feof(influence_report_file)) {
                scanf_error (cmodel, i, attr_idx, influence_report_filename);
                return NULL;
            }
        }

        cmodel->class_weight[i] /= n_cases;

        /* Looks for the beginning of the description list */
        while ((fscanf(influence_report_file,
                      "%*s REAL ATTRIBUTE (t = R) |Mean-jk - "
                       "numb t mtt description "
                       "I-jk Mean StDev Mean-*k|/ Mean StDev "
                       "t a -jk -jk StDev-jk -*k -*%c",
                       &dummy_c) != 1) || (dummy_c != 'k')) {
            if (feof(influence_report_file)) {
                scanf_error (cmodel, i, attr_idx, influence_report_filename);
                return NULL;
            }
        }
        i_ = cmodel->dim;
        while (i_){
            if (fscanf (influence_report_file, "%u %u %*s %s %5c",
                        (unsigned int*)&corr_idx, (unsigned int*)&j,
                        dummy_str, dummy_s ) != 4){
               scanf_error (cmodel, i, attr_idx, influence_report_filename);
               return NULL;
            }
            if (dummy_s[0] == 'p'){
                if (fscanf (influence_report_file,
                            "%*f ( %e %e ) %*e ( %*e %*e )",
                            &mean, &dev ) != 2){
                    scanf_error (cmodel, i, attr_idx,
                                 influence_report_filename);
                    return NULL;
                }
                /* Maps the attribute index to the storage index */
                j = (size_t *)bsearch(&j, attr_idx, n_attr,sizeof(attr_idx[0]),
                                      compare) -
                    (size_t *)attr_idx;
                cmodel->param[i][j][0] = mean; /* Save mean and deviation */
                cmodel->param[i][j][1] = dev;

                /*  decides which model this dimension uses */
                if (strcmp(dummy_str, "MNcn") == 0) {
                    cmodel->classmodel[i][j] = MULTI_NORMAL;
                    cmodel->param[i][j][2] = corr_idx;
                }
                else if (strcmp(dummy_str, "SNcn") == 0)
                    cmodel->classmodel[i][j] = SINGLE_NORMAL;
                else
                    cmodel->classmodel[i][j] = UNKNOWN;
                i_--;
            }
            else if (dummy_s[0] =='a'){
                fscanf (influence_report_file, "%*f %*s %*e %*e %*e");
                for (k = 1; k < MIN_AA; k++) {
                    if (fscanf( influence_report_file, "%s %*e %*e %*e",
                                dummy_str) !=1){
                        scanf_error (cmodel, i, attr_idx,
                                     influence_report_filename);
                        return NULL;
                    }
                }
            }
            else {
                scanf_error (cmodel, i, attr_idx, influence_report_filename);
                return NULL;
            }
        }

        /* Tries to read the correlation matrices */
        fgetpos (influence_report_file, &pos);
        fscanf (influence_report_file, "%s", dummy_str);
        while (strcmp("DATA_CORR_MATRIX", dummy_str)) {
            if (feof (influence_report_file)) {
                scanf_error (cmodel, i, attr_idx, influence_report_filename);
                return NULL;
            }
            fgetpos (influence_report_file, &pos);
            fscanf (influence_report_file, "%s", dummy_str);
        }
        /* set the pos to the beginning of the "DATA_CORR_MATRIX" */
        fsetpos (influence_report_file, &pos);
        while (fscanf(influence_report_file, "%s ", dummy_str) == 1
               && strcmp(dummy_str, "DATA_CORR_MATRIX") == 0) {
            fgetpos(influence_report_file, &pos);

            while (fscanf(influence_report_file, "%*s %u",
                          (unsigned int*)&dummy_i) != 1) {
                fgetpos(influence_report_file, &pos);
                pos__ = pos_ = pos;
            }

            /* maps the attribute index to the storage index */
            j = (size_t*)bsearch(&dummy_i, attr_idx,
                                 n_attr, sizeof(attr_idx[0]), compare)
                - (size_t*)attr_idx;
            use_idx[j] = 1;
            j__ = j_ = dummy_i;

            while (fscanf(influence_report_file, "%u",
                          (unsigned int*)&dummy_i) == 1) {
                j__ = j_;     /* saves previously read numbers to have the */
                j_ = dummy_i; /* right ones by hand */
                /* maps the attribute index to the storage index */
                j = (size_t*)bsearch(&j__, attr_idx,
                                     n_attr, sizeof(attr_idx[0]), compare)
                    - (size_t*)attr_idx;
                use_idx[j] = 1;
                pos__ = pos_;            /* saves previous file positions to */
                pos_ = pos;              /* to jump back later */
                fgetpos(influence_report_file, &pos);
            }
            /* Maps the attribute index to the storage index */
            j = (size_t*)bsearch(&j__, attr_idx,
                                 n_attr, sizeof(attr_idx[0]), compare)
                - (size_t*)attr_idx;
            use_idx[j] = 1;
            pos = pos__;

            /* Reads the correlation values for the found attributes */
            fsetpos(influence_report_file, &pos);
            for (k = 0; k < cmodel->dim; k++) {
                if (use_idx[k] == 0)
                    continue;
                fscanf(influence_report_file, "%u", (unsigned int*)&dummy_i);

                for (l = 0; l < cmodel->dim; l++) {
                    if (use_idx[l] == 0)
                        continue;
                    fscanf(influence_report_file, "%lf",
                           &cmodel->cov_matrix[i][(k*cmodel->dim)+l]);
                    cmodel->cov_matrix[i][(k*cmodel->dim)+l] *=
                        cmodel->param[i][k][1]*cmodel->param[i][l][1];
                }
            }
            for (k = 0; k < cmodel->dim; k++)
                use_idx[k] = 0;
        }
        fsetpos(influence_report_file, &pos);
    }

    free_if_not_null(use_idx);
    free_if_not_null(attr_idx);
    fclose(influence_report_file);

    return cmodel;
}


/* ---------------- computeMemebership -------------------------
 * For a case descriptor test_vec (of dimension cmodel->dim)
 * compute the memberships to each class in cmodel.
 * Returns the vector of memberships (length is cmodel->n_class)
 * or NULL if an unsupported class is found in cmodel.
 */
/* 27 MAR 2006 TStehr
 * This function was broken into two parts:
 * computeMembership and computeMemberStrct
 * It shouldn't affect any routine which calls computeMembership
 */
/* --------------- computeMembershipStrct ----------------------
 * Core of the original function computeMembershipStrct
 * mship should be initialized, normally to class_weight or 1
 */
float *
computeMembershipStrct (float *mship, const float * test_vec,
                        const struct clssfcn *cmodel)
{
    const char this_sub[] = "computeMembershipStrct";
    size_t i, j , k;

    if (mship == NULL) {
        err_printf (this_sub, "mship is null\n");
        return NULL;
    }
    if (test_vec == NULL) {
       err_printf (this_sub, "test_vec is NULL\n");
       return NULL;
    }
		
    for (i = 0; i < cmodel->n_class; i++) {   /*for every class...*/
        switch (cmodel->classmodel[i][0]) {
        case SINGLE_NORMAL:
			mprintf("single normal\n");
            for (j = 0; j < cmodel->dim; j++) {
                mship[i] *= (1/cmodel->param[i][j][1])*
                    (gsl_sf_erf_Q((test_vec[j] - fabs(cmodel->abs_error)
                                   - cmodel->param[i][j][0])
                                  / cmodel->param[i][j][1]) -
                     gsl_sf_erf_Q((test_vec[j] + fabs(cmodel->abs_error)
                                   - cmodel->param[i][j][0])
                                  / cmodel->param[i][j][1]));
            }
            break;
        case MULTI_NORMAL: {
			mprintf("multi normal\n");
            double *C_matrix;
            gsl_matrix_view C;
            gsl_permutation *p;
            gsl_matrix* C_inv;
            double C_det;
            double *y_array;
            gsl_vector_view y;
            gsl_vector* C_inv_y;
            double exponent = 1.0;
            size_t to_mall;
            int s;
            
            to_mall = cmodel->dim * cmodel->dim;
            C_matrix = E_CALLOC(to_mall, sizeof (C_matrix[0]));
            for (k = 0; k < to_mall; k++)
                C_matrix[k] = cmodel->cov_matrix[i][k];
            
            
            C = gsl_matrix_view_array(C_matrix, cmodel->dim, cmodel->dim);
            p = gsl_permutation_alloc(cmodel->dim);
            
            gsl_linalg_LU_decomp(&C.matrix, p, &s);
            
            C_inv = gsl_matrix_alloc(cmodel->dim, cmodel->dim);
            gsl_linalg_LU_invert(&C.matrix, p, C_inv);
            
            C_det = gsl_linalg_LU_det(&C.matrix, s);
            
            y_array = E_CALLOC(cmodel->dim, sizeof(double));
            for (k = 0; k < cmodel->dim; k++)
                if (test_vec[k] == BAD_ANGLE)
                    y_array[k] = 0.0;
                else
                    y_array[k] = test_vec[k] - cmodel->param[i][k][0];
            
            y = gsl_vector_view_array(y_array,cmodel->dim);
            
            C_inv_y = gsl_vector_alloc(cmodel->dim);
            
            gsl_blas_dgemv(CblasNoTrans, 1.0, C_inv, &y.vector, 0.0, C_inv_y);
            
            gsl_blas_ddot(C_inv_y, &y.vector, &exponent);
            
            exponent *= -0.5;
            mship[i] *= 2 * cmodel->abs_error *
                exp(exponent) / sqrt(/* 2 * M_PI *  */fabs(C_det));
            
            gsl_vector_free(C_inv_y);
            gsl_matrix_free(C_inv);
            free_if_not_null(y_array);
            free_if_not_null(C_matrix);
            gsl_permutation_free(p);
        }
            break;
        case UNKNOWN:
        default:
            err_printf(this_sub, "Unknown classmodel\n");
            return NULL;
        }
    }
    return mship;
}

float *
computeMembership(float *mship, const float* test_vec,
                  const struct clssfcn *cmodel)
{
    long double sum = 0.0;
    size_t i;
    
    for (i = 0; i < cmodel->n_class; i++) {
        mship[i] = cmodel->class_weight[i];   /* Default class probability */
    }
    if (computeMembershipStrct (mship, test_vec, cmodel) == NULL)
        return NULL;
    
    for (i = 0; i < cmodel->n_class; i++)
        sum += (mship[i] * mship[i]);
    sum = sqrt (sum);
    if (sum != 0)
        for (i = 0; i < cmodel->n_class; i++)
            mship[i] /= sum;
    return mship;
}

/***************************************************************
float *
computeMembership(float *mship, const float* test_vec,
                  const struct clssfcn *cmodel)
{
    const char this_sub[] = "computeMembership";
    long double sum = 0.0;
    size_t i, j, k;

    if (mship == NULL) {
        err_printf(this_sub, "mship is null\n");
        return NULL;
    }

    if (test_vec == NULL) {
        err_printf(this_sub, "test_vec is NULL\n");
        return NULL;
    }

    for (i = 0; i < cmodel->n_class; i++) {
        mship[i] = cmodel->class_weight[i];
        switch (cmodel->classmodel[i][0]) {
        case SINGLE_NORMAL:
            for (j = 0; j < cmodel->dim; j++) {
                mship[i] *= (1/cmodel->param[i][j][1])*
                    (gsl_sf_erf_Q((test_vec[j] - fabs(cmodel->abs_error)
                                   - cmodel->param[i][j][0])
                                  / cmodel->param[i][j][1])
                     - gsl_sf_erf_Q((test_vec[j] + fabs(cmodel->abs_error)
                                     - cmodel->param[i][j][0])
                                    / cmodel->param[i][j][1]));
            }
            break;
        case MULTI_NORMAL: {
            double *C_matrix;
            gsl_matrix_view C;
            gsl_permutation *p;
            gsl_matrix* C_inv;
            double C_det;
            double *y_array;
            gsl_vector_view y;
            gsl_vector* C_inv_y;
            double exponent = 1.0;
            size_t to_mall;
            int s;

            to_mall = cmodel->dim * cmodel->dim;
            C_matrix = E_CALLOC(to_mall, sizeof (C_matrix[0]));
            for (k = 0; k < cmodel->dim*cmodel->dim; k++)
                C_matrix[k] = cmodel->cov_matrix[i][k];


            C = gsl_matrix_view_array(C_matrix, cmodel->dim, cmodel->dim);
            p = gsl_permutation_alloc(cmodel->dim);

            gsl_linalg_LU_decomp(&C.matrix, p, &s);

            C_inv = gsl_matrix_alloc(cmodel->dim, cmodel->dim);
            gsl_linalg_LU_invert(&C.matrix, p, C_inv);

            C_det = gsl_linalg_LU_det(&C.matrix, s);

            y_array = E_CALLOC(cmodel->dim, sizeof(double));
            for (k = 0; k < cmodel->dim; k++)
                if (test_vec[k] == BAD_ANGLE)
                    y_array[k] = 0.0;
                else
                    y_array[k] = test_vec[k] - cmodel->param[i][k][0];

            y = gsl_vector_view_array(y_array,cmodel->dim);

            C_inv_y = gsl_vector_alloc(cmodel->dim);
            gsl_blas_dgemv(CblasNoTrans, 1.0, C_inv, &y.vector, 0.0, C_inv_y);

            gsl_blas_ddot(C_inv_y, &y.vector, &exponent);

            exponent *= -0.5;
            mship[i] *=
                2*cmodel->abs_error * (1/sqrt(fabs(C_det)))*exp(exponent);

            gsl_vector_free(C_inv_y);
            gsl_matrix_free(C_inv);
            free_if_not_null(y_array);
            free_if_not_null(C_matrix);
            gsl_permutation_free(p);
        }
            break;
        case UNKNOWN:
        default:
            err_printf(this_sub, "Unknown classmodel\n");
            return NULL;
        }
        sum += (mship[i] * mship[i]);
    }
    sum = sqrt (sum);
    if (sum != 0.0)
        for (i = 0; i < cmodel->n_class; i++)
            mship[i] /= sum;


    return mship;
}
***************************************************************/
