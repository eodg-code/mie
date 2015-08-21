;+
 PRO mie_lognormal, Nd, Rm, Sg, Wavenumber, Cm, Dqv=dqv, dlm=dlm, $
                    npts=npts, xres=xres, info=info, $;mthread=mthread, $
                    Bext, Bsca, w, g, ph, Bbac=Bbac

; NAME:
;       MIE_LOGNORMAL
; PURPOSE:
;       Calculates the scattering parameters of a log normal
;       distribution of spherical particles.
;
; CATEGORY:
;       EODG Mie routines
; CALLING SEQUENCE:
;       mie_lognormal, Nd, Rm, Sg, Wavenumber, Cm $
;       ,Bext ,Bsca ,w ,g ,ph [, Dqv = dqv] [, Dlm = Dlm] $
;       [, Bbac = Bbac] 
; INPUTS:
;       Nd:         Number density of the particle distribution
;       Rm:         Median radius of the particle distribution
;       Sg:         The spread of the distribution, such that the
;                   standard deviation of ln(r) is ln(S)
;       Wavenumber: Wavenumber of light (units must match Rm)
;       Cm:         Complex refractive index
; OPTIONAL KEYWORD PARAMETERS:
;       Dqv:        An array of the cosines of scattering angles at
;                   which to compute the phase function.
;       dlm:        If set the IDL DLM version of the Mie scattering
;                   procedure will be called rather than the IDL coded
;                   version.
;       npts:       Allows the user to overide the automatically
;                   calculated number of quadrature points for the
;                   integration over size. NOTE: reducing the number of
;                   abscissa can substantially decrease the accuracy of
;                   the result - BE CAREFUL!!!!
;       xres:       Sets the spacing of the quadrature points (in size
;                   parameter). Overridden by npts. Default is 0.1. The
;                   same warning as npts applies here!
;       info:       Named variable that, on return, will contain a
;                   structure containing the number of abscissa points
;                   and the maximum and minimum size parameters used.
;       mthread:    Controls the number of threads which will be
;                   utilised by the DLM version Mie code.
;                   HAS NO EFFECT UNLESS dlm IS ALSO SET
;                   By default the code will use up to 4 threads. The
;                   behaviour of the code for different values of this
;                   keyword is as follows:
;                   mthread=0 : Code will utilise only 1 thread.
;                   mthread=1 : Code will utilise the number of threads
;                               specified by the !CPU.TPOOL_NTHREADS
;                               IDL system variable.
;                   mthread>1 : Code wiil utilise the number of threads
;                               specifed by mthread, or the total
;                               number available on the system
;                               (whichever is less).
;                   * BE AWARE * : Simply adding more threads will not
;                     necessarily speed up the calculation beyond a
;                     factor of about 4, and may decrease the accuracy
;                     of the result! It's recommended that no more
;                     than the default number of threads are used
;                     without testing!                           
;
; OUTPUT PARAMETERS:
;       Bext:       The extinction coefficient
;       Bsca:       The scattering coefficient
;       w:          The single scatter albedo
;       g:          The asymmetry parameter
;       ph:         The phase function - an array of the same
;                   dimension as Dqv. Only calculated if Dqv is
;                   specified.
;       Bbac:       The backscatter coefficient
; RESTRICTIONS:
;       Note, this procedure calls the MIE_SINGLE (or MIE_DLM_SINGLE),
;       QUADRATURE and SHIFT_QUADRATURE procedures.
; MODIFICATION HISTORY
;       G. Thomas Sept. 2003 (based on mie.pro written by Don Grainger)
;
;       G. Thomas Nov. 2003 minor bug fixes
;
;       G. Thomas Feb. 2004 Explicit double precission added throughout
;       G. Thomas Apr. 2005 DLM keyword added to enable the use of mie_dlm_single
;       RGG       Jun  2005 Implemented 0.1 step size in X and trapezium quadrature
;       G. Thomas Jun. 2005 Added calculation of the phase function afer
;       calling mie_dlm_single, as the DLM no longer returns it. Changed "size" to
;       "Dx" (as size is a IDL keyword!). Also added npts and info keywords.
;       RGG        8 Jun 2005 Slight modification of code.
;       G. Thomas  9 Jun 2005 Added xres keyword
;       A. Smith   6 Jul 2009 Added backscatter coefficient and +ve cm warning.
;       G. Thomas 21 Jul 2011 Added mthread keyword, to make use of the
;                             parallelised version of the Mie DLM code
;       G. Thomas 22 Jul 2011 Can't get parallelised mie DLM to work
;                             correctly; S1 and S2 arrays corrupted in
;                             the Fortran. Giving up - mthread keyword
;                             code commented out
;       G. Thomas 12 Jun 2012 "Bug" (more an IDL language "feature" really) fix:
;                             If the number density was passed as a single
;                             element array, size distribution was being
;                             truncated to smallest size.
;-
    Common mieln, absc, wght

    ru_max = 10000d0
;   Check the Nd only has one element
    if n_elements(Nd) gt 1 then message,'Number density should be a scalar quantity!'

;   Create vectors for size integration
    Tq = gauss_cvf(0.999D0)
    
    Rl = exp(alog(Rm)+Tq*alog(Sg))
    Ru = exp(alog(Rm)-Tq*alog(Sg)+alog(4))

    if 2D0 * !dpi * Rl * Wavenumber ge ru_max then message,'Lower bound of integral is larger than maximum permitted size parameter.'
    if 2D0 * !dpi * Ru * Wavenumber ge ru_max then begin
        Ru = (ru_max - 1d0) / ( 2D0 * !dpi * Wavenumber )
        message,/continue,'Warning: Radius upper bound truncated to avoid size parameter overflow.'
    endif

  if imaginary(cm) gt 0d0 and not(keyword_set(silent)) then $
    message, /continue,'Warning: Imaginary part of refractive index '+$
             'should be negative for absorbing particles. '+$
             'Set /silent to hide this message.'

    if not keyword_set(xres) then xres = 0.1

    if not keyword_set(npts) then begin
;   Accurate calulation requires 0.1 step size in x

       Npts = (Long(2D0 * !dpi * (ru-rl) * Wavenumber/xres)) > 200
    endif

;   quadrature on the radii
    if n_elements(wght) ne Npts then quadrature,'T',Npts,absc,wght

    shift_quadrature,absc,wght,Rl,Ru,R,W1

    W1P = Nd[0] * R * W1 * sqrt(!dpi)* $
          EXP(-0.5D0*(ALOG(R/Rm) / ALOG(Sg))^2) / (sqrt(2D0) * ALOG(Sg))

    Dx = 2D0 * !dpi * R * Wavenumber

    if arg_present(info) then info = { Npts    : Npts, $
                                       MinSize : Dx[0], $
                                       MaxSize : Dx[Npts-1] }

;   If an array of cos(theta) is provided, calculate phase function
    if n_elements(dqv) gt 0 then begin
        Inp = n_elements(dqv)
        ph = dblarr(Inp)
        if keyword_set(dlm) then begin
;          Put the mthread keyword into the right form for the DLM call...
;           if n_elements(mthread) gt 0 then begin
;              if mthread eq 0 then mthrd = 0 $
;              else if mthread eq 1 then mthrd = !CPU.TPOOL_NTHREADS $
;              else if mthread gt 1 then mthrd = mthread $
;              else mthrd = 4
;           endif else mthrd = 4
           DCm = dcomplex(Cm)   ; Ensure the inputs are double precision
           DDqv = double(Dqv)
           Mie_dlm_single, Dx, DCm, Dqv=DDqv, Dqxt, Dqsc, Dqbk, Dg, $
                           Xs1, Xs2;, mthread=mthrd
;          Cannot get the DLM to return the phase function. Very misterious...
           Dph = dblarr(Inp,Npts)
           for i = 0,Inp-1 do $
              Dph[i,*] = 2d0 * double(Xs1[i,*]*CONJ(Xs1[i,*]) + Xs2[i,*]*CONJ(Xs2[i,*])) $
              / (Dx^2 * Dqsc)
        endif else begin
;          Put the mthread keyword into the right form for the DLM call...
;           if n_elements(mthread) gt 0 then begin
;              if mthread eq 0 then mthrd = 0 $
;              else if mthread eq 1 then mthrd = !CPU.TPOOL_NTHREADS $
;              else if mthread gt 1 then mthrd = mthread $
;              else mthrd = 4
;           endif else mthrd = 4
           DCm = dcomplex(Cm)   ; Ensure the inputs are double precision
           DDqv = double(Dqv)
           Mie_single, Dx, DCm, Dqv=DDqv, Dqxt, Dqsc, Dqbk, Dg, $
                       Xs1, Xs2, Dph;, mthread=mthrd
        endelse
    endif else begin
        inp = 1
        dqv = 1D0
        if keyword_set(dlm) then begin
            Mie_dlm_single, double(Dx),dcomplex(Cm),Dqxt,Dqsc,Dqbk,Dg
        endif else begin
            Mie_single, Dx,Cm,Dqxt,Dqsc,Dqbk,Dg
        endelse
    endelse

    Bext = total(W1P * DQxt)
    Bsca = total(W1P * DQsc)
    g = total(W1P * dg * DQsc) / Bsca
    w = Bsca / Bext
    Bbac = total(W1P * DQbk)

    if n_elements(ph) gt 0 then $ ;Phase function
      for i =0,Inp-1 do $
        ph(i) = total(W1P * Dph(i,*) * Dqsc) / Bsca
END
