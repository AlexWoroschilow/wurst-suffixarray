/*
 * 26 March 2002
 * rscid = $Id: vec_i.h,v 1.1 2002/04/02 00:04:16 torda Exp $
 */

#ifndef VEC_I_H
#define VEC_I_H

struct RPoint *
vector_difference (struct RPoint *res, const struct RPoint *a,
                   const struct RPoint *b);
struct RPoint *
vector_add (struct RPoint *res, const struct RPoint *a,
            const struct RPoint *b);

float scalar_product(const struct RPoint *a, const struct RPoint *b);

struct RPoint *
vec_scl (struct RPoint *res, const struct RPoint *v, const float a);

float vector_length (const struct RPoint *v);
float vector_sqr_length (const struct RPoint *v);

struct RPoint *
vec_nrm (struct RPoint *res, const struct RPoint *v, const float r);

struct RPoint *
vector_product (struct RPoint *res, const struct RPoint *u,
                const struct RPoint *v);

#endif
