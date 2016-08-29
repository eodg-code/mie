;+
; NAME:
;     mie_derivs_ln
;
; PURPOSE:
;     Calculates the scattering parameters and their analytical
;     deritatives (wrt the parameters of the distribution) of a log
;     normal distribution of spherical particles.
;
;     A the derivation of expressions for the analytical derivatives of
;     Mie scattering terms is covered by:
;     Grainger, R.G., J. Lucas, G.E. Thomas, G. Ewan, "The Calculation
;     of Mie Derivatives", Submitted to Appl. Opt., 2004.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;     mie_derivs_ln, N, Rm, S, Wavenumber, Cm[, Dqv=dqv] $
;     [, Bext][, Bsca][, dBextdN][, dBextdRm][, dBextdS] $
;     [, dBscadN][, dBscadRm][, dBscadS][, i1][, i2] $
;     [, di1dN][, di1dRm][, di1dS][, di2dN][, di2dRm][, di2dS]
;
; INPUTS:
;     N:          The number density of the particle distribution
;     Rm:         Median radius of the particle distribution (microns)
;     S:          The spread of the distribution, such that the stardard
;                 deviation of ln(r) is ln(S)
;     Wavenumber: Wavenumber of light, defined as 1/wavelength. A
;                 positive scalar whos units match those of Rm.
;     Cm:         Complex refractive index
;
; KEYWORD PARAMETERS:
;     Dqv:        An array of the cosines of scattering angles at which
;                 to compute the intensity functions.
;
; OUTPUTS:
;     Bext:       The volume extinction coefficient
;     Bsca:       The volume scattering coefficient
;     dBextdN:    Derivative of the extinction coefficient wrt the
;                 number density
;     dBextdRm:   Derivative of the extinction coefficient wrt the
;                 mean radius
;     dBextdS:    Derivative of the extinction coefficient wrt the
;                 spread
;     dBscadN:    Derivative of the scattering coefficient wrt the
;                 number density
;     dBscadRm:   Derivative of the scattering coefficient wrt the
;                 mean radius
;     dBscadS:    Derivative of the scattering coefficient wrt the
;                 spread
;     i1:         The first intensity function - intensity of light
;                 polarised in the plane perpendicular to the directions
;                 of inicident light propogation and observation. Only
;                 caculated if dqv is specified.
;     i2:         The second intensity function - intensity of light
;                 polarised in the plane parallel to the directions of
;                 inicident light propogation and observation. Only
;                 caculated if dqv is specified.
;     di1dN:      Derivatives of the first intensity function wrt the
;                 number density
;     di1dRm:     Derivatives of the first intensity function wrt the
;                 mean radius
;     di1dS:      Derivatives of the first intensity function wrt the
;                 spread
;     di2dN:      Derivatives of the second intensity function wrt the
;                 number density
;     di2dRm:     Derivatives of the second intensity function wrt the
;                 mean radius
;     di2dS:      Derivatives of the second intensity function wrt the
;                 spread
;                 NB. The values of involving the intensity functions
;                 are arrays of the same dimension as dqv, and are only
;                 calculated if dqv is specified.
;
; KEYWORD OUTPUTS:
;
; RESTRICTIONS:
;     Note, this procedure calls the MIE_SINGLE and QUADRATURE
;     procedures.
;
; MODIFICATION HISTORY:
;     G. Thomas, Sep 2003: mie_derivs_ln.pro
;     G. Thomas, Nov 2003: minor bug fixes
;     G. Thomas, Feb 2004: Explicit double precission added throughout
;         and header updated.
;     G. Thomas, Jun 2005: Implemented 0.1 step size in X and trapezium
;         quadrature. Also added npts and info keywords.
;-

pro mie_derivs_ln, N,Rm,S,Wavenumber,Cm,Dqv=dqv,Bext,Bsca, $
                   dBextdN,dBextdRm,dBextdS,dBscadN,dBscadRm,dBscadS, $
                   i1,i2,di1dN,di1dRm,di1dS,di2dN,di2dRm,di2dS

    Common miedervln, absc, wght

;   Create vectors for size integration
    Tq = gauss_cvf(0.999D0)
    Rl = exp(alog(Rm)+Tq*alog(S))  & Ru = exp(alog(Rm)-Tq*alog(S)+alog(4))
    if 2D0 * !dpi * Ru * Wavenumber ge 12000 then begin
        Ru = 11999D0 / ( 2D0 * !dpi * Wavenumber )
        print,'mie_lognormal: Warning! Radius upper bound truncated to avoid size parameter overflow.'
    endif

    if not keyword_set(npts) then begin
;      Accurate calulation requires 0.1 step size in x
       Npts = (Long(2D0 * !dpi * (ru-rl) * Wavenumber/0.1)) > 200
    endif

;   quadrature on the radii
    if n_elements(wght) ne Npts then quadrature,'T',Npts,absc,wght

    shift_quadrature,absc,wght,Rl,Ru,R,Wghtr

    Tmp =  EXP(-0.5D0*(ALOG(R/Rm) / ALOG(S))^2) / (sqrt(2D0*!dpi) * ALOG(S) * R)

    W1 = N * Tmp

    Dx = 2D0 * !dpi * R * Wavenumber
    if keyword_set(info) then info = { Npts    : Npts, $
                                       MinSize : Dx[0], $
                                       MaxSize : Dx[Npts-1] }

;   Create Mie variables
    Dqxt=DBLARR(npts)
    Dqsc=DBLARR(npts)

    Dx = 2d0 * !dpi * R * Wavenumber
;   If an array of cos(theta) is provided, calculate phase function
    if n_elements(dqv) gt 0 then begin
        Inp = n_elements(dqv)
        i1 = dblarr(Inp) & i2 = i1
        di1dN = i1 &  di2dN = i1
        di1dRm = i1 & di2dRm = i1
        di1dS = i1 & di2dS = i1
        Mie_single, Dx,Cm,Dqv=dqv,Dqxt,Dqsc,Dqbk,Dg,Xs1,Xs2,Dph
    endif else begin
        inp = 1
        dqv = 1d0
        Mie_single, Dx,Cm,Dqxt,Dqsc,Dqbk,Dg
    endelse

    lnRRm = ALOG(R/Rm) ;Precalculate for speed

    Bext = Total(wghtr * W1 * DQxt * !dpi * R^2) ; Extinction
    dBextdN = Bext / N
    dBextdRm = Total(wghtr * W1 * DQxt * lnRRm * !dpi * R^2) / (alog(S)^2 * Rm)
    dBextdS = (Total(wghtr * W1 * DQxt * lnRRm^2 * !dpi * R^2) / alog(S)^2 $
             - Total(wghtr * W1 * DQxt * !dpi * R^2)) / (S * alog(S))

    Bsca = Total(wghtr * W1 * DQsc * !dpi * R^2) ;Scattering
    dBscadN = Total(wghtr * W1 * DQsc * !dpi * R^2) / N
    dBscadRm = Total(wghtr * W1 * DQsc * lnRRm * !dpi * R^2) / (alog(S)^2 * Rm)
    dBscadS = Total(wghtr * W1 * DQsc * lnRRm^2 * !dpi * R^2) / (S * alog(S)^3) $
            - Total(wghtr * W1 * DQsc * !dpi * R^2) / (S * alog(S))

    if n_elements(i1) gt 0 then $ ; Intensity functions
        for i =0,Inp-1 do begin
            i1(i) = Total(wghtr * W1 * real_part(Xs1(i,*)*conj(Xs1(i,*))))
            di1dN(i) = i1(i) / N
            di1dRm(i) = Total(wghtr * W1 * lnRRm * real_part(Xs1(i,*)*conj(Xs1(i,*)))) $
                        / (alog(S)^2 * Rm)
            di1dS(i) = Total(wghtr * W1 * lnRRm * real_part(Xs1(i,*)*conj(Xs1(i,*)))) $
                       / (S * alog(S)^3) $
                       - Total(wghtr * W1 * real_part(Xs1(i,*)*conj(Xs1(i,*)))) $
                       / (S * alog(S))

            i2(i) = Total(wghtr * W1 * real_part(Xs2(i,*)*conj(Xs2(i,*))))
            di2dN(i) = i2(i) / N
            di2dRm(i) = Total(wghtr * W1 * lnRRm * real_part(Xs2(i,*)*conj(Xs2(i,*)))) $
                        / (alog(S)^2 * Rm)
            di2dS(i) = Total(wghtr * W1 * lnRRm * real_part(Xs2(i,*)*conj(Xs2(i,*)))) $
                       / (S * alog(S)^3) $
                       - Total(wghtr * W1 * real_part(Xs2(i,*)*conj(Xs2(i,*)))) $
                       / (S * alog(S))
    endfor
end
