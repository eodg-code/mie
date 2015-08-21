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

/* ANSI */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <omp.h>

/* IDL */
#include "export.h"

/* Local */
#include "mie_dlm_single-openmp.h"

/* FOR IDL >= 5.3 */
/* Define the procedures */
static IDL_SYSFUN_DEF2 mie_single_procedures[] = {
  {(IDL_FUN_RET) mie_dlm_single,"MIE_DLM_SINGLE",2,10,IDL_SYSFUN_DEF_F_KEYWORDS,0},
};

/* Startup call when DLM is loaded */
int IDL_Load(void)
{
    /* Add procedures */
    if (!IDL_SysRtnAdd(mie_single_procedures, FALSE,
                    ARRLEN(mie_single_procedures))) {
        return IDL_FALSE;
    }


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
    IDL_LONG64 dims[2], tdims[2];
    int tid;
    int tthread, nthread;
    /* Variables to be passed to and from Mie routine itself */
    IDL_LONG mie_ok;
    IDL_LONG Npts;
    IDL_LONG Inp;
    IDL_DCOMPLEX Cm;
    IDL_DCOMPLEX  *S1ARR, *S2ARR;
    double Qext, Qsca, Qbsc, g;
    double *S1r, *S1i, *S2r, *S2i, *Phase;
    double *DqvARR;
    /* Variables to hold values for different sizes */
    double *DxARR, *QextARR, *QscaARR, *QbscARR, *gARR, *PhaseARR;
    /* Arrays for passing back to IDL */
    IDL_VPTR ivQextArr, ivQscaArr, ivQbscArr, ivgArr, ivS1, ivS2, ivPhase,
             ivS1Arr, ivS2Arr, ivPhaseArr;
    /* Definition of keyword parameters Dqv  - description stored in Dqvdesc */
    /*                                      - data stored in Dqvdata */
    typedef struct {
      IDL_KW_RESULT_FIRST_FIELD;
      int mthread;
      int DqvExists;
      IDL_MEMINT DqvOffset;
      double DqvData[10000];
    } KW_RESULT;
    
    static IDL_KW_ARR_DESC_R DqvDesc = {IDL_KW_OFFSETOF(DqvData), 1, 10000, 
					IDL_KW_OFFSETOF(DqvOffset)};
    static IDL_KW_PAR kw_pars[] = {
      IDL_KW_FAST_SCAN,
      {"DQV",IDL_TYP_DOUBLE,1,IDL_KW_ARRAY,IDL_KW_OFFSETOF(DqvExists),
       IDL_CHARA(DqvDesc)},
      {"MTHREAD",IDL_TYP_INT,1,IDL_KW_ZERO,0,IDL_KW_OFFSETOF(mthread)},
      {NULL}};
    KW_RESULT kw;
    printf("argc :%d\n",argc);
    argc = IDL_KWProcessByOffset(argc, argv, argk, kw_pars, (IDL_VPTR *)0, 
				 1, &kw);
    printf("argc :%d\n",argc);
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

    /* Set up parallel processing of the "mieint" loop:
       How many threads to use. If the mthread keyword has been set
       to one, or is greater than the total number of avaialble 
       threads, use all threads offered by the computer. If it's value
       is zero or negative, only use one thread, otherwise, use its
       value */
    tthread = omp_get_max_threads();
    if (kw.mthread > tthread || kw.mthread == 1) nthread = tthread;
    else if (kw.mthread < 1) nthread = 1; else nthread = kw.mthread;
    omp_set_num_threads(nthread);
    /* printf("%d threads available, %d threads used, based on a mthread value of %d\n",
       tthread, nthread, kw.mthread); */

    if (kw.DqvExists)
      {
	Inp = kw.DqvOffset;
	dims[0] = Inp; dims[1] = Npts;
	/* Dqv is defined, so use the user supplied scattering angles */
	DqvARR = (double *)malloc(sizeof(double[Inp]));
	for (j=0; j<Inp; j++) DqvARR[j] = kw.DqvData[j];

	tdims[0] = Inp; tdims[1] = nthread;
	/* Set up the output array temporaries */
	S1ARR    = (IDL_DCOMPLEX *) IDL_MakeTempArray(IDL_TYP_DCOMPLEX, 2, dims,
						      IDL_ARR_INI_ZERO, &ivS1Arr);
	S2ARR    = (IDL_DCOMPLEX *) IDL_MakeTempArray(IDL_TYP_DCOMPLEX, 2, dims,
						      IDL_ARR_INI_ZERO, &ivS2Arr);
	PhaseARR = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE, 2, dims,
						IDL_ARR_INI_ZERO, &ivPhaseArr);
	
	/* Initalise the parallel region with a compiler construct to the openmp library */
#pragma omp parallel num_threads(nthread)			\
  shared(DxARR, Cm, Inp, DqvARR, QextARR, QscaARR, QbscARR,	\
	 gARR, S1ARR, S2ARR, PhaseARR)				\
  private(i, j, mie_ok, Qext, Qsca, Qbsc, g, S1r, S1i,		\
	  S2r, S2i, Phase, tid)
	{	  
	  /* Set up temporary storage arrays for the output from mieint. Note that
	     is a bit messy - as far as I can see the IDL API doesn't provide a
	     way of cleaning up array temporaries made using it's own constructors.
	     Therefore, I generate arrays which a explicitly indexed by thread
	     number (rather than just setting them as private variables */
	  S1r   = (double *)malloc(sizeof(double[Inp]));
	  S1i   = (double *)malloc(sizeof(double[Inp]));
	  S2r   = (double *)malloc(sizeof(double[Inp]));
	  S2i   = (double *)malloc(sizeof(double[Inp]));
	  Phase = (double *)malloc(sizeof(double[Inp]));
	  
#pragma omp for schedule(guided)
	  for ( i = 0; i < Npts; i++)
	    {
	      tid = omp_get_thread_num();
	      /* Call the mieint procedure */
	      mieintnocmplx_(&DxARR[i], &Cm, &Inp, DqvARR, &Qext, 
			     &Qsca, &Qbsc, &g, S1r, S1i, S2r, S2i, 
			     Phase, &mie_ok);
	      if (mie_ok != 0)
		IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
			    "Size parameter overflow!");

	      printf("S1: ");
	      for ( j = 0; j < Inp; j++ ) printf("(%f, %f) ",S1r[i],S1i[i]);
	      printf("\n");

	      QextARR[i] = Qext;
	      QscaARR[i] = Qsca;
	      QbscARR[i] = Qbsc;
	      gARR[i] = g;
	      for ( j = 0; j < Inp; j++ )
		{
		  S1ARR[i*Inp+j].r = S1r[j];
		  S1ARR[i*Inp+j].i = S1i[j];
		  S2ARR[i*Inp+j].r = S2r[j];
		  S2ARR[i*Inp+j].i = S2i[j];
		  PhaseARR[i*Inp+j] = Phase[j];
		}
	      
	    }
	printf("\n");
	} /* End of Parallel code section */
	
        /* Check to see if each of the output parameters has been
           specified in the IDL call, if it has copy the value to
           the return argument, otherwise simply deallocate it.
           NOTE: The keyword is included as an argument in defining
           argc but NOT argv - crazy!*/
	
        if (argc > 8)
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
    else
      {
	/* Dqv has not been defined, so set Inp = 1 and do calculations
	   at 0 degrees */
	Inp = 1;
	DqvARR = (double *)malloc(sizeof(double));
	DqvARR[0] = 1;
	tdims[0] = Inp; tdims[1] = nthread;
	/* Set up the output array temporaries */
	S1ARR    = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX, Npts,
						       IDL_ARR_INI_ZERO, &ivS1Arr);
	S2ARR    = (IDL_DCOMPLEX *) IDL_MakeTempVector(IDL_TYP_DCOMPLEX, Npts,
						       IDL_ARR_INI_ZERO, &ivS2Arr);
	PhaseARR = (double *) IDL_MakeTempVector(IDL_TYP_DOUBLE, Npts, 
						 IDL_ARR_INI_ZERO, &ivPhaseArr);
	
	/* Initalise the parallel region with a compiler construct to the openmp library */
#pragma omp parallel num_threads(nthread)			\
  shared(DxARR, Cm, Inp, DqvARR, QextARR, QscaARR, QbscARR,	\
	 gARR, S1ARR, S2ARR, PhaseARR)				\
  private(i, j, mie_ok, Qext, Qsca, Qbsc, g, S1r, S1i, S2r, S2i,\
	  Phase, tid)
	{
	  /* Set up temporary storage arrays for the output from mieint. Note that
	     is a bit messy - as far as I can see the IDL API doesn't provide a
	     way of cleaning up array temporaries made using it's own constructors.
	     Therefore, I generate arrays which a explicitly indexed by thread
	     number (rather than just setting them as private variables */
	  S1r   = (double *)malloc(sizeof(double[Inp]));
	  S1i   = (double *)malloc(sizeof(double[Inp]));
	  S2r   = (double *)malloc(sizeof(double[Inp]));
	  S2i   = (double *)malloc(sizeof(double[Inp]));
	  Phase = (double *)malloc(sizeof(double[Inp]));

#pragma omp for schedule(guided)
	  for ( i = 0; i < Npts; i++)
	    {
	      tid = omp_get_thread_num();
	      /* Call the mieint procedure */
	      mieintnocmplx_(&DxARR[i], &Cm, &Inp, DqvARR, &Qext, 
			     &Qsca, &Qbsc, &g, S1r, S1i, S2r, S2i, 
			     Phase, &mie_ok);
	      if (mie_ok != 0)
		IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP,
			    "Size parameter overflow!");
	      
	      QextARR[i] = Qext;
	      QscaARR[i] = Qsca;
	      QbscARR[i] = Qbsc;
	      gARR[i] = g;
	      S1ARR[i].r = S1r[0];
	      S1ARR[i].i = S1i[0];
	      S2ARR[i].r = S2r[0];
	      S2ARR[i].i = S2i[0];
	      PhaseARR[i] = Phase[0];
	    }
	} /* End of Parallel code section */
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
    
    /* Reset the max number of threads to the original values */
    omp_set_num_threads(tthread);
    
    /* Cleanup any temporaries due to the keyword and
       deallocate the remain temporary arrays.*/
    IDL_KW_FREE;
    
    free(S1r);    S1r = NULL;
    free(S1i);    S1i = NULL;
    free(S2r);    S2r = NULL;
    free(S2i);    S2i = NULL;
    free(DqvARR); DqvARR = NULL;
    
    return;
}
