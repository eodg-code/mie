; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


;+
; NAME:
;     mie_size_dist
;
; PURPOSE:
;     Calculates the scattering parameters of a size distribution of spherical
;     particles.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;     mie_size_dist, distname, Nd, params, wavenumber, Cm [, Dqv = Dqv] $
;     [, Npts=Npts] [, xres=xres] [, /DLM] [, mthread=mthread] [, /SILENT], $
;     Bext, Bsca, w, g [, SPM] [, Bbac=Bbac] [, Gavg=Gavg] [, Vavg=Vavg] $
;     [, Ravg=Ravg] [, RVW=RVW] [, info=info]
;
; INPUTS:
;     distname    Name of the size distribution. 'log_normal' or 'modified_gamma'
;                 according to Hansen and Travis 1974.
;     Nd:         Number density of the particle distribution
;     params:     Array of size distribution parameters
;                 distname eq 'log_normal'
;                     params[0] : Median radius of the particle distribution
;                     params[1] : The spread of the distribution, such that the
;                                 standard deviation of ln(r) is ln(S)
;                 distname eq 'modified_gamma'
;                     params[0] : a
;                     params[1] : b
;                     params[2] : minimum radius in the distribution
;                     params[3] : maximum radius in the distribution
;     wavenumber: Wavenumber of light (units must match units of the size
;                 distribution)
;     Cm:         Complex refractive index of the particle(s)
;
; KEYWORD PARAMETERS:
;     Dqv:        Cosines of scattering angles at which to compute the intensity
;                 functions etc.
;     Npts:       If set, this keyword overrides the default calculation of the
;                 quadrature points (which provides points at 0.1 spacing in the
;                 size parameter). NOTE: reducing the number of abscissa can
;                 substantially decrease the accuracy of the result. BE CAREFUL!
;     xres:       Sets the spacing of the quadrature points (in size parameter).
;                 Overridden by Npts. Default is 0.1. The same warning as Npts
;                 applies here!
;     DLM:        If set the IDL DLM version of the algorithm will be used
;                 instead of the IDL coded version.
;     mthread:    Controls the number of threads which will be utilised by the
;                 DLM version of the algorithm. If not set by default the code
;                 will use 1 thread. The behaviour of the code for different
;                 values of this keyword is as follows:
;                 mthread<=0 : !CPU.TPOOL_NTHREADS
;                 mthread> 0 : Code will utilise the number of threads specified
;                              by mthread.
;                 * BE AWARE * : Adding more threads than the number of physical
;                   cores (hyperthreads do not count as physical cores) on the
;                   system will not speed up the calculation.
;
;                 * HAS NO EFFECT UNLESS DLM IS ALSO SET*
;     SILENT:     If set all warning messages issued by the code will be
;                 suppressed.
;
; OUTPUTS:
;     Bext:       Total extinction coefficient of the distribution
;     Bsca:       Total scattering coefficient of the distribution
;     w:          Single scatter albedo of the distribution
;     g:          Asymmetry parameter of the distribution
;
; OPTIONAL OUTPUTS:
;     SPM:        Scattering phase matrix elements F11 (SPM[0,*]), F33
;                 (SPM[1,*]), F12 (SPM[2,*]), F34 (SPM[3,*]), where the 2nd
;                 dimension is the same dimension as Dqv. Also only calculated
;                 if Dqv is specified.
;
; KEYWORD OUTPUTS:
;     Bbac:       Backscatter coefficient
;     Gavg:       Average projected area per particle
;     Vavg:       Average volume per particle
;     Ravg:       Average radius
;     RVW:        Volume-weighted average radius
;     info:       Named variable that, on return, will contain a structure
;                 containing the number of abscissa points and the maximum and
;                 minimum size parameters used.
;
; RESTRICTIONS:
;     Note, this procedure calls the mie_single (or mie_dlm_single), quadrature
;     and shift_quadrature procedures.
;
; MODIFICATION HISTORY:
;     G. Thomas, Sep 2003: Based on mie.pro written by Don Grainger.
;     G. Thomas, Nov 2003: Minor bug fixes
;     G. Thomas, Feb 2004: Explicit double precision added throughout.
;     G. Thomas, Apr 2005: DLM keyword added to enable the use of
;         mie_dlm_single.
;     R. Grainger, Jun 2005: Implemented 0.1 step size in X and trapezium
;         quadrature
;     G. Thomas, Jun 2005: Added calculation of the phase function after calling
;         mie_dlm_single, as the DLM no longer returns it. Changed "size" to
;         "Dx" (as size is a IDL keyword!). Also added Npts and info keywords.
;     R. Grainger, 8 Jun 2005: Slight modification of code.
;     G. Thomas, 9 Jun 2005: Added xres keyword
;     A. Smith,  6 Jul 2009: Added backscatter coefficient and +ve cm warning.
;     G. Thomas, 21 Jul 2011: Added mthread keyword, to make use of the
;         parallelised version of the Mie DLM code
;     G. Thomas, 22 Jul 2011: Can't get parallelised mie DLM to work correctly;
;         S1 and S2 arrays corrupted in the Fortran. Giving up - mthread keyword
;         code commented out
;     G. Thomas, 12 Jun 2012: "Bug" (more an IDL language "feature" really) fix:
;         If the number density was passed as a single element array, size
;         distribution was being truncated to smallest size.
;     G. McGarragh, 29 Jul 2015: Added initial support for additional size
;         distributions starting with the gamma distribution as defined in
;         Hansen and Travis 1974.
;     G. McGarragh, 29 Jul 2015: Changed scattering phase function output to
;         scattering phase matrix output.
;     G. McGarragh, 29 Jul 2015: DLM output of the phase matrix SPM was fixed so
;         no need to calculate it here any more.
;     G. McGarragh, 29 Jul 2015: Added support to optionally output several
;         geometric parameters averaged over the size distribution including
;         Gavg, Vavg, Ravg, and RVW.
;     G. McGarragh, 10 Dec 2015: Better naming of size distributions: gamma ->
;         modified_gamma and lognormal -> log_normal.
;-

pro mie_size_dist, distname, Nd, params, wavenumber, Cm, Dqv=Dqv, Npts=Npts, $
                   xres=xres, info=info, DLM=DLM, mthread=mthread, $
                   SILENT=SILENT, Bext, Bsca, w, g, SPM, Bbac=Bbac, Gavg=Gavg, $
                   Vavg=Vavg, Ravg=Ravg, RVW=RVW

    Common mieln, absc, wght

    ru_max = 10000d0
;   Check the Nd only has one element
    if n_elements(Nd) gt 1 then message,'Number density should be a scalar quantity!'

;   Create vectors for size integration
    Tq = gauss_cvf(0.999D0)

    if distname eq 'modified_gamma' then begin
        Rl = params[2]
        Ru = params[3]
    endif else if distname eq 'log_normal' then begin
        Rl = exp(alog(params[0])+Tq*alog(params[1]))
        Ru = exp(alog(params[0])-Tq*alog(params[1])+alog(4))
    endif else begin
        message,'Invalid size distribution name: ' + distname
    endelse

    if 2D0 * !dpi * Rl * wavenumber ge ru_max then message,'Lower bound of ' + $
        'integral is larger than maximum permitted size parameter.'
    if 2D0 * !dpi * Ru * wavenumber ge ru_max then begin
        Ru = (ru_max - 1d0) / ( 2D0 * !dpi * wavenumber )
        message,/continue,'Warning: Radius upper bound truncated to avoid '  + $
            'size parameter overflow.'
    endif

    if imaginary(cm) gt 0d0 and not(keyword_set(SILENT)) then $
        message, /continue,'Warning: Imaginary part of refractive index '+ $
            'should be negative for absorbing particles. Set /SILENT to '+ $
            'hide this message.'

    if not keyword_set(xres) then xres = 0.1

    if not keyword_set(Npts) then begin
;       Accurate calculation requires 0.1 step size in x
        Npts = (long(2D0 * !dpi * (ru-rl) * wavenumber/xres)) > 200
    endif

;   quadrature on the radii
    if n_elements(wght) ne Npts then quadrature,'T',Npts,absc,wght

    shift_quadrature,absc,wght,Rl,Ru,R,W1

    if distname eq 'modified_gamma' then begin
        W1P = W1 * Nd[0] * R^((1. - 3. * params[1]) / params[1]) * $
              exp(-R / (params[0] * params[1]))
    endif else if distname eq 'log_normal' then begin
        W1P = W1 * Nd[0] / (sqrt(2D0) * sqrt(!dpi) * R * alog(params[1])) * $
              exp(-0.5D0*(alog(R/params[0]) / alog(params[1]))^2)
    endif

    Dx = 2D0 * !dpi * R * wavenumber

    if arg_present(info) then info = { Npts    : Npts, $
                                       MinSize : Dx[0], $
                                       MaxSize : Dx[Npts-1] }

;   If an array of cos(theta) is provided, calculate phase matrix
    if n_elements(Dqv) gt 0 then begin
        Inp = n_elements(Dqv)
        SPM = dblarr(4,Inp)
        if keyword_set(DLM) then begin
;           Put the mthread keyword into the right form for the DLM call...
            if keyword_set(mthread) gt 0 then begin
                if mthread lt 1 then mthrd = !CPU.TPOOL_NTHREADS $
                else mthrd = mthread
            endif else mthrd = 1
            DCm = dcomplex(Cm) ; Ensure the inputs are double precision
            DDqv = double(Dqv)
            Mie_dlm_single, Dx, DCm, Dqv=DDqv, Dqxt, Dqsc, Dqbk, Dg, Xs1, Xs2, $
                            F11, F33, F12, F34 ;, mthread=mthrd
            DSPM = dblarr(4,Inp,Npts)
            DSPM[0,*,*] = F11
            DSPM[1,*,*] = F33
            DSPM[2,*,*] = F12
            DSPM[3,*,*] = F34
;           Mie_dlm_single, Dx, DCm, Dqv=DDqv, Dqxt, Dqsc, Dqbk, Dg, Xs1, Xs2 ; $
;                           ;, mthread=mthrd
;           Cannot get the DLM to return the phase matrix. Very misterious...
;           So must calculate the elements of DSPM below.
;           DSPM = dblarr(4,Inp,Npts)
;           AA = 2d0 / (Dx^2 * Dqsc)
;           for i = 0,Inp-1 do begin
;               DSPM[0,i,*] =  AA * double( Xs1[i,*]*conj(Xs1[i,*]) + $
;                                           Xs2[i,*]*conj(Xs2[i,*]))
;               DSPM[1,i,*] =  AA * double( Xs1[i,*]*conj(Xs2[i,*]) + $
;                                           Xs2[i,*]*conj(Xs1[i,*]))
;               DSPM[2,i,*] = -AA * double( Xs1[i,*]*conj(Xs1[i,*]) - $
;                                           Xs2[i,*]*conj(Xs2[i,*]))
;               DSPM[3,i,*] = -AA * double((Xs1[i,*]*conj(Xs2[i,*]) - $
;                                           Xs2[i,*]*conj(Xs1[i,*])) * $
;                                   complex(0.0d, 1.0d))
;           endfor
        endif else begin
;           Put the mthread keyword into the right form for the DLM call...
            if keyword_set(mthread) gt 0 then begin
                if mthread lt 1 then mthrd = !CPU.TPOOL_NTHREADS $
                else mthrd = mthread
            endif else mthrd = 1
            DCm = dcomplex(Cm) ; Ensure the inputs are double precision
            DDqv = double(Dqv)
            Mie_single, Dx, DCm, Dqv=DDqv, Dqxt, Dqsc, Dqbk, Dg, Xs1, Xs2, $
                        DSPM ;, mthread=mthrd
        endelse
    endif else begin
        if keyword_set(DLM) then begin
            Mie_dlm_single, double(Dx), dcomplex(Cm), Dqxt, Dqsc, Dqbk, Dg
        endif else begin
            Mie_single, Dx, Cm, Dqxt, Dqsc, Dqbk, Dg
        endelse
    endelse

    W1PA = W1P * !dpi * R^2
    W1PV = W1PA * 4./3. * R

    Bext = total(W1PA * DQxt)
    Bsca = total(W1PA * DQsc)
    g = total(W1PA * dg * DQsc) / Bsca
    w = Bsca / Bext
    if arg_present(Bbac) then Bbac = total(W1PA * DQbk)

    if n_elements(Dqv) gt 0 then begin
        for i = 0,Inp-1 do begin
            SPM[0,i] = total(W1PA * DSPM[0,i,*] * Dqsc) / Bsca
            SPM[1,i] = total(W1PA * DSPM[1,i,*] * Dqsc) / Bsca
            SPM[2,i] = total(W1PA * DSPM[2,i,*] * Dqsc) / Bsca
            SPM[3,i] = total(W1PA * DSPM[3,i,*] * Dqsc) / Bsca
        endfor
    endif

    if arg_present(Gavg) then Gavg = total(W1PA)
    if arg_present(Vavg) then Vavg = total(W1PV)
    if arg_present(Ravg) then Ravg = total(W1P * R)
    if arg_present(RVW) then RVW = total(W1PV * R) / Vavg
end
