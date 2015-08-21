/* File:      mie_dlm_single.c
 * Purpose:   Wrapper for the FORTRAN routine mieint
              in the file mieint.f
 * Author:    Gareth Thomas
 * Version:   2.0
 * Date:      July 2011
 * Known bugs: Returning the phase function has been disabled.
               For some reason if 9 variables are returned to IDL
               a segmentation fault results. I've got no idea why 
* Hisotry:    Version 1.0 June 2005.
              Version 2.0 July 2011.
	                - Added the mthread keyword parameter and
			  corresponding openmp compiler directives
			  and functions.
			- Changed to the new(er) IDL_KWProcessByOffset
			  function for dealing with Keywords (as old
			  version isn't supported by export.h in
			  IDL >8.0). Also slightly changed the way
			  the Dqv array is handled.
                        - Removed support for IDL before version
                          5.3, as new keyword handler wouldn't have
                          worked anyway.
*/

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
extern void mieintnocmplx_(const double *Dx, const IDL_DCOMPLEX *Cm, const IDL_LONG *inp, const double *Dqv, double *Qext, double *Qsca, double *g, double *Qbsc, double *S1r, double *S1i, double *S2r, double *S2i, double *Phase,  IDL_LONG *mie_ok);
#endif
