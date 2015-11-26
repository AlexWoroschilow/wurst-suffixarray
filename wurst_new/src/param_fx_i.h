/*
 * 10 Jan 2002
 * rcsid = $Id: param_fx_i.h,v 1.3 2002/04/09 07:29:48 torda Exp $
 */

#ifndef PARAM_FX_I_H
#define PARAM_FX_I_H
struct FXParam;
struct FXParam *param_fx_read (const char *fname);
void            FXParam_destroy (struct FXParam *fx);
#endif  /* PARAM_FX_I_H */
