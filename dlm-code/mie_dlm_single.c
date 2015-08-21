/* File:      mie_dlm_single.c
 * Purpose:   Wrapper for the FORTRAN routine mieint
              in the file mieint.f
 * Author:    Gareth Thomas
 * Version:   1.1
 * Date:      June 2005
 * Known bugs: Returning the phase function has been disabled.
               For some reason if 9 variables are returned to IDL
               a segmentation fault results. I've got no idea why */

/* ANSI */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

/* IDL */
#include "export.h"

/* Local */
#include "mie_dlm_single.h"

#ifdef __IDLPRE53__
    /* FOR IDL < 5.3 */
    /* Define the procedures */
    static IDL_SYSFUN_DEF mie_single_procedures[] = {
        {(IDL_FUN_RET) mie_dlm_single,"MIE_DLM_SINGLE",2,10,IDL_SYSFUN_DEF_F_KEYWORDS},
    };
#else
    /* FOR IDL >= 5.3 */
    /* Define the procedures */
    static IDL_SYSFUN_DEF2 mie_single_procedures[] = {
        {(IDL_FUN_RET) mie_dlm_single,"MIE_DLM_SINGLE",2,10,IDL_SYSFUN_DEF_F_KEYWORDS,0},
    };
#endif


/* Startup call when DLM is loaded */
int IDL_Load(void)
{
    /* IDL version 5.3 and greater use IDL_SYSFUN_DEF2 while earlier versions use
     * IDL_SYSFUN_DEF.  Note the addition of the final '0' in each line for
     * IDL_SYSFUN_DEF2. */

    /* If IDL is pre-5.3 then change IDL_SysRtnAdd to IDL_AddSystemRoutine,
     * (NB: the parameters stay the same) */

#ifdef __IDLPRE53__
    /* FOR IDL < 5.3 */
    /* Add procedures */
    if (!IDL_AddSystemRoutine(mie_single_procedures, FALSE,
                    ARRLEN(mie_single_procedures))) {
                return IDL_FALSE;
        }

#else
    /* FOR IDL >= 5.3 */
    /* Add procedures */
    if (!IDL_SysRtnAdd(mie_single_procedures, FALSE,
                    ARRLEN(mie_single_procedures))) {
        return IDL_FALSE;
    }

#endif

    /* Register the error handler */
    IDL_ExitRegister(mie_single_exit_handler);
    return(IDL_TRUE);
}

/* Called when IDL is shutdown */
void mie_single_exit_handler(void)
{
/* Nothing special to do in this case */
}

/* Wrapped Fortran routines below this point */

/* =======================================================================
 * IDL Procedure:  mie_dlm_single
 * Description:    Performs Mie computation for an array of particle sizes
 *                 for a given refractive index. Expected values are:
 *                 Dx:        Array of size parameters (MUST be an ARRAY of
 *                            type DOUBLE, even if it has one element!!!)
 *                 Cm:        Complex refractive index (MUST be be a scalar
 *                            of type DOUBLE COMPLEX)
 *                 Keywords:
 *                 Dqv = Dqv  An ARRAY of cos(scattering angles) at which
 *                            intensity functions and phase function are
 *                            evaluated. If not specified, a single value
 *                            of 1.0 (forward scatter) is used.
 *                 Returns:
 *                 Qext:      Array of extinction efficiencies
 *                 Qsca:      Array of scattering efficiencies
 *                 Qbsc:      Array of back-scatter efficiencies
 *                 g:         Array of asymetry parameters
 *                 s1 and s2: Arrays of scattering functions (in theta
 *                            and size)
 *                 Phase:     Array of Phase function (one for each size)
 * G Thomas November 2003
 * =======================================================================
 */
void IDL_CDECL mie_dlm_single(int argc, IDL_VPTR argv[], char *argk)
{
    /* local */
    long i, j;
    long dims[2];
    double dqvtmp[1];
    /* Variables to be passed to and from Mie routine itself */
    IDL_LONG mie_ok;
    IDL_LONG Npts;
    IDL_LONG Inp;
    IDL_DCOMPLEX Cm;
    IDL_DCOMPLEX  *S1, *S2, *S1ARR, *S2ARR;
    double Qext, Qsca, Qbsc, g;
    double *Phase;
    /* Variables to hold values for different sizes */
    double *DxARR, *QextARR, *QscaARR, *QbscARR, *gARR, *PhaseARR;
    /* Arrays for passing back to IDL */
    IDL_VPTR ivQextArr, ivQscaArr, ivQbscArr, ivgArr, ivS1, ivS2, ivPhase,
             ivS1Arr, ivS2Arr, ivPhaseArr;
    /* Definition of keyword parameter Dqv  - description stored in Dqvdesc */
    /*                                      - data stored in Dqvdata */
    IDL_VPTR outargv[1];
    static double Dqvdata[10000];
    static int dqvexists;
    static IDL_KW_ARR_DESC Dqvdesc = {(char*)Dqvdata,1,10000,0};
    static IDL_KW_PAR kw_pars[] = {
                IDL_KW_FAST_SCAN,
                {"DQV",IDL_TYP_DOUBLE,1,IDL_KW_ARRAY,&dqvexists,IDL_CHARA(Dqvdesc)},
                {NULL}};
    /* Ensure all temporary stuff resulting from keyword is cleaned up */
    IDL_KWCleanup(IDL_KW_MARK);
    /* Extraction of keyword and normal parameters */
    IDL_KWGetParams(argc,argv,argk,kw_pars,outargv,1);
    /* Ensure the refractive index is given as a scalar */
    IDL_ENSURE_SCALAR(argv[1]);

    if (argv[0]->type != IDL_TYP_DOUBLE)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
          "Size parameters must be DOUBLE precision");
    if (argv[1]->type != IDL_TYP_DCOMPLEX)
    IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
          "Refractive index must be DOUBLE precision COMPLEX");
    IDL_ENSURE_ARRAY(argv[0]);

    /* assign values to local variables */
    mie_ok   = 0;
    Npts     = argv[0]->value.arr->n_elts;
    DxARR    = (double *) argv[0]->value.arr->data;
    Cm       = argv[1]->value.dcmp;
    QextARR  = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                               IDL_ARR_INI_ZERO, &ivQextArr);
    QscaARR  = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                                 IDL_ARR_INI_ZERO, &ivQscaArr);
    QbscARR  = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                                 IDL_ARR_INI_ZERO, &ivQbscArr);
    gARR     = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts,
                                 IDL_ARR_INI_ZERO, &ivgArr);

    /* Copy the complex refractive index value to the C
         complex<double> variable Cm_c */

  if (dqvexists)
    {
        /* Dqv is defined, so use the user supplied scattering angles */
        Inp = Dqvdesc.n;
        dims[0] = Inp; dims[1] = Npts;
        S1       = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                     Inp, IDL_ARR_INI_ZERO, &ivS1);
        S2       = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                     Inp, IDL_ARR_INI_ZERO, &ivS2);
        Phase    = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
                                     Inp, IDL_ARR_INI_ZERO, &ivPhase);
        S1ARR    = (IDL_DCOMPLEX *) IDL_MakeTempArray(IDL_TYP_DCOMPLEX, 2,
                                     dims, IDL_ARR_INI_ZERO, &ivS1Arr);
        S2ARR    = (IDL_DCOMPLEX *) IDL_MakeTempArray(IDL_TYP_DCOMPLEX, 2,
                                     dims, IDL_ARR_INI_ZERO, &ivS2Arr);
        PhaseARR = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, 2,
                                     dims, IDL_ARR_INI_ZERO, &ivPhaseArr);

        for ( i = 0; i < Npts; i++)
        {
            /* Call the mieint procedure */
            mieint_(&DxARR[i],&Cm,&Inp,Dqvdata,&Qext,&Qsca,&Qbsc,&g,S1,S2,Phase,&mie_ok);
            if (mie_ok != 0)
               IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                           "Size parameter overflow!");
            QextARR[i] = Qext;
            QscaARR[i] = Qsca;
            QbscARR[i] = Qbsc;
            gARR[i] = g;
            for ( j = 0; j < Inp; j++ )
            {
                S1ARR[i*Inp+j] = S1[j];
                S2ARR[i*Inp+j] = S2[j];
                PhaseARR[i*Inp+j] = Phase[j];
            }
        }
        /* Check to see if each of the output parameters has been
           specified in the IDL call, if it has copy the value to
           the return argument, otherwise simply deallocate it.
           NOTE: The keyword is included as an argument in defining
           argc buy NOT argv - crazy!*/

        if (argc > 9)
        {
           /*IDL_VarCopy(ivPhaseArr, argv[8]);*/
           IDL_VarCopy(ivS2Arr, argv[7]);
           IDL_VarCopy(ivS1Arr, argv[6]);
           IDL_VarCopy(ivgArr, argv[5]);
           IDL_VarCopy(ivQbscArr, argv[4]);
           IDL_VarCopy(ivQscaArr, argv[3]);
           IDL_VarCopy(ivQextArr,  argv[2]);
        }
        else
        {
           IDL_DELTMP(ivPhaseArr);
           if (argc > 8)
           {
              IDL_VarCopy(ivS2Arr, argv[7]);
              IDL_VarCopy(ivS1Arr, argv[6]);
              IDL_VarCopy(ivgArr, argv[5]);
              IDL_VarCopy(ivQbscArr, argv[4]);
              IDL_VarCopy(ivQscaArr, argv[3]);
              IDL_VarCopy(ivQextArr,  argv[2]);
           }
           else
           {
              IDL_DELTMP(ivS2Arr);
              if (argc > 7)
              {
                 IDL_VarCopy(ivS1Arr, argv[6]);
                 IDL_VarCopy(ivgArr, argv[5]);
                 IDL_VarCopy(ivQbscArr, argv[4]);
                 IDL_VarCopy(ivQscaArr, argv[3]);
                 IDL_VarCopy(ivQextArr,  argv[2]);
              }
              else
              {
                 IDL_DELTMP(ivS1Arr);
                 if (argc > 6)
                 {
                    IDL_VarCopy(ivgArr, argv[5]);
                    IDL_VarCopy(ivQbscArr, argv[4]);
                    IDL_VarCopy(ivQscaArr, argv[3]);
                    IDL_VarCopy(ivQextArr,  argv[2]);
                 }
                 else
                 {
                    IDL_DELTMP(ivgArr);
                    if (argc > 5)
                    {
                       IDL_VarCopy(ivQbscArr, argv[4]);
                       IDL_VarCopy(ivQscaArr, argv[3]);
                       IDL_VarCopy(ivQextArr,  argv[2]);
                    }
                    else
                    {
                       IDL_DELTMP(ivQbscArr);
                       if (argc > 4)
                       {
                          IDL_VarCopy(ivQscaArr,  argv[3]);
                          IDL_VarCopy(ivQextArr,  argv[2]);
                       }
                       else
                       {
                          IDL_DELTMP(ivQscaArr);
                          if (argc > 3)
                             IDL_VarCopy(ivQextArr,  argv[2]);
                          else
                             IDL_DELTMP(ivQextArr)
                       }
                    }
                 }
              }
           }
        }
    }
    else
    {
        /* Dqv has not been defined, so set Inp = 1 and do calculations
           at 0 degrees */
        Inp = 1;
        dqvtmp[0] = 1;
        dims[0] = Inp; dims[1] = Npts;
        S1       = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                     Inp, IDL_ARR_INI_ZERO, &ivS1);
        S2       = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                     Inp, IDL_ARR_INI_ZERO, &ivS2);
        Phase    = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
                                     Inp, IDL_ARR_INI_ZERO, &ivPhase);
        S1ARR    = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                     Npts, IDL_ARR_INI_ZERO, &ivS1Arr);
        S2ARR    = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX,
                                     Npts, IDL_ARR_INI_ZERO, &ivS2Arr);
        PhaseARR = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE,
                                     Npts, IDL_ARR_INI_ZERO, &ivPhaseArr);

        for ( i = 0; i < Npts; i++)
        {
            /* Call the mieint procedure */
            mieint_(&DxARR[i],&Cm,&Inp,dqvtmp,&Qext,&Qsca,&Qbsc,&g,S1,S2,Phase,&mie_ok);
            if (mie_ok != 0)
               IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
                           "Size parameter overflow!");

            QextARR[i] = Qext;
            QscaARR[i] = Qsca;
            QbscARR[i] = Qbsc;
            gARR[i] = g;
            S1ARR[i] = S1[0];
            S2ARR[i] = S2[0];
            PhaseARR[i] = Phase[0];
        }
        /* Check to see if each of the output parameters has been
           specified in the IDL call, if it has copy the value to
           the return argument, otherwise simply deallocate it. */

        if (argc > 8)
        {
           IDL_VarCopy(ivPhaseArr, argv[8]);
           IDL_VarCopy(ivS2Arr, argv[7]);
           IDL_VarCopy(ivS1Arr, argv[6]);
           IDL_VarCopy(ivgArr, argv[5]);
           IDL_VarCopy(ivQbscArr, argv[4]);
           IDL_VarCopy(ivQscaArr, argv[3]);
           IDL_VarCopy(ivQextArr,  argv[2]);
        }
        else
        {
           IDL_DELTMP(ivPhaseArr);
           if (argc > 7)
           {
              IDL_VarCopy(ivS2Arr, argv[7]);
              IDL_VarCopy(ivS1Arr, argv[6]);
              IDL_VarCopy(ivgArr, argv[5]);
              IDL_VarCopy(ivQbscArr, argv[4]);
              IDL_VarCopy(ivQscaArr, argv[3]);
              IDL_VarCopy(ivQextArr,  argv[2]);
           }
           else
           {
              IDL_DELTMP(ivS2Arr);
              if (argc > 6)
              {
                 IDL_VarCopy(ivS1Arr, argv[6]);
                 IDL_VarCopy(ivgArr, argv[5]);
                 IDL_VarCopy(ivQbscArr, argv[4]);
                 IDL_VarCopy(ivQscaArr, argv[3]);
                 IDL_VarCopy(ivQextArr,  argv[2]);
              }
              else
              {
                 IDL_DELTMP(ivS1Arr);
                 if (argc > 5)
                 {
                    IDL_VarCopy(ivgArr, argv[5]);
                    IDL_VarCopy(ivQbscArr, argv[4]);
                    IDL_VarCopy(ivQscaArr, argv[3]);
                    IDL_VarCopy(ivQextArr,  argv[2]);
                 }
                 else
                 {
                    IDL_DELTMP(ivgArr);
                    if (argc > 4)
                    {
                       IDL_VarCopy(ivQbscArr, argv[4]);
                       IDL_VarCopy(ivQscaArr, argv[3]);
                       IDL_VarCopy(ivQextArr,  argv[2]);
                    }
                    else
                    {
                       IDL_DELTMP(ivQbscArr);
                       if (argc > 3)
                       {
                          IDL_VarCopy(ivQscaArr,  argv[3]);
                          IDL_VarCopy(ivQextArr,  argv[2]);
                       }
                       else
                       {
                          IDL_DELTMP(ivQscaArr);
                          if (argc > 2)
                             IDL_VarCopy(ivQextArr,  argv[2]);
                          else
                             IDL_DELTMP(ivQextArr)
                       }
                    }
                 }
              }
           }
        }
    }

  /* Cleanup any temporaries due to the keyword and
     deallocate the remain temporary arrays.*/
  IDL_KWCleanup(IDL_KW_CLEAN);

  IDL_DELTMP(ivS1);
  IDL_DELTMP(ivS2);
  IDL_DELTMP(ivPhase);
  return;
}
