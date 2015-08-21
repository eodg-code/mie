/* File:      mie_dlm_single.h
 * Purpose:   Header for mie_dlm_single.c
 * Author:    Gareth Thomas
 * Version:   1.1
 * Date:      June 2005
 * Known bugs: Returning the phase function has been disabled.
               For some reason if 9 variables are returned to IDL
               a segmentation fault results. I've got no idea why */
#ifndef IDL_C_miesingle
#define IDL_C_miesingle

/* Message Numbers */
#define mie_single_ERROR         0
#define mie_single_NOSTRINGARRAY    -1

/* Useful macro */
#define ARRLEN(arr) (sizeof(arr)/sizeof(arr[0]))

extern IDL_MSG_BLOCK msg_block;

/* Define the startup function that adds C functions to IDL
 * along with the exit handler */

/* IDL-Postres interface */
extern void mie_single_exit_handler(void);
extern int  mie_single_Startup(void);

/* Define the wrapper function which calls the fortran function */
extern void IDL_CDECL mie_dlm_single(int argc, IDL_VPTR argv[], char *argk);
extern void mieint_(const double *Dx, const IDL_DCOMPLEX *Cm, const IDL_LONG *inp, const double *Dqv, double *Qext, double *Qsca, double *g, double *Qbsc, IDL_DCOMPLEX *S1, IDL_DCOMPLEX *S2, double *Phase,  IDL_LONG *mie_ok);
#endif
