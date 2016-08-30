;+
; NAME:
;     mie_derivs
;
; PURPOSE:
;     Calculates the scattering parameters and their analytical derivatives of a
;     series of particles using Mie scattering theory.
;
;     A the derivation of expressions for the analytical derivatives of Mie
;     scattering terms is covered by:
;     Grainger, R.G., J. Lucas, G.E. Thomas, G. Ewan, "The Calculation of Mie
;     Derivatives", Submitted to Appl. Opt., 2004.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;     mie_derivs, x, Cm, Dqv $
;     [, Qext] [, Qsca] [, dQextdx] [, dQextdRem] [, dQextdImm] $
;     [, dQscadx] [, dQscadRem] [, dQscadImm] [, i1] [, i2] $
;     [, di1dx] [, di2dx] [, di1dRem] [, di1dImm] [, di2dRem] [, di2dImm] $
;     [, ASYM=asym] [, dASYMdx=dasymdx] [, dASYMdRem=dasymdrem] $
;     [, dASYMdImm=dASYMdrem] [, /SILENT]
;
; INPUTS:
;     x:         The particle size parameter
;     Cm:        The complex refractive index of the particle
;     Dqv:       The cosine of the scattering angles at which to calculate the
;                intensity functions etc
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;     SILENT:    Don't give warnings about positive k values.
;
; OUTPUTS:
;     Qext:      The extinction efficiency
;     Qsca:      The extinction efficiency
;     dQextdx:   Derivative of the extinction efficiency wrt the particle size
;                parameter
;     dQextdRem: Derivative of the extinction efficiency wrt the real part of
;                the refractive index
;     dQextdImm: Derivative of the extinction efficiency wrt the imaginary part
;                of the refractive index
;     dQscadx:   Derivative of the scattering efficiency wrt the particle size
;                parameter
;     dQscadRem: Derivative of the scattering efficiency wrt the real part of
;                the refractive index
;     dQscadImm: Derivative of the scattering efficiency wrt the imaginary part
;                of the refractive index
;     i1:        The first intensity function - intensity of light polarized in
;                the plane perpendicular to the directions of incident light
;                propagation and observation.
;     i2:        The second intensity function - intensity of light polarized in
;                the plane parallel to the directions of incident light
;                propagation and observation.
;     di1dx:     Derivatives of the intensity functions wrt the
;     di2dx:     particle size parameter
;     di1dRem:   Derivatives of the first intensity function wrt to
;     di1dImm:   the real and imaginary parts of the refractive index
;     di2dRem:   Derivatives of the second intensity function wrt
;     di2dImm:   to the real and imaginary parts of the refractive index
;
;     NB. The values of involving the intensity functions are arrays of the same
;     dimension as dqv and are only calculated if dqv is specified.
;
; OPTIONAL OUTPUTS:
;
; KEYWORD OUTPUTS:
;     asym:      The asymmetry parameter.
;     dasymdx:   Derivative of asymmetry wrt size parameter.
;     dasymdRem: Derivative of asymmetry wrt real pt of RI.
;     dasymdImm: Derivative of asymmetry wrt imag pt of RI.
;
; RESTRICTIONS:
;
; MODIFICATION HISTORY:
;     J. Lucas, 2002: Main line programme
;     R. Grainger, 2002: Procedure version
;     G. Thomas, Sep 2003: Put into EODG routines format
;     J. Graham (UC Berkeley), Feb 2004: Introduced explicit double precision
;         numerical values into all computational expressions.
;     G. Thomas, Feb 2004: Header updated.
;     A. Smith, Aug 2010: Added positive k warning.
;     A. Smith, Apr 2013: Added asymmetry parameter.
;-

pro mie_derivs, x,Cm,Dqv,Qext,Qsca, $
                dQextdx,dQextdRem,dQextdImm, $
                dQscadx,dQscadRem,dQscadImm, $
                i1,i2,di1dx,di2dx,di1dRem,di1dImm,di2dRem,di2dImm, $
                silent=silent, $
                asym=asym,dasymdx=dasymdx, $
                dasymdRem=dasymdRem,dasymdImm=dasymdImm

    if imaginary(cm) gt 0d0 and not(keyword_set(silent)) then $
        message,/continue,'Warning: Imaginary part of refractive index '+ $
            'should be negative for absorbing particles. Set /silent to '+ $
            'hide this message.'

    if x lt 0.02 then NStop = 2 else begin
        case 1 OF
            x le 8.0    : NStop = x + 4.00*x^(1./3.) + 1.0
            x lt 4200.0 : NStop = x + 4.05*x^(1./3.) + 2.0
            ; no. of terms required for convergence of Mie expressions
            ; (Wiscombe 1980)
            else        : NStop = x + 4.00*x^(1./3.) + 2.0
        endcase
    endelse

    m = cm

    y = m*x
    Nmx = fix(max([NStop,abs(y)]) + 15.)
    D = dcomplexarr(Nmx+1)

    for p = Nmx-1,1,-1 do begin
        R = (p+1) / y
        D(p) = R - 1/(R+D(p+1)) ; downward recurrence to find A_n(y)
    endfor

    psim1 = cos(x)
    psi0 = sin(x)
    chim1 = -sin(x)
    chi0 = cos(x)

    Pi0 = 0.D0
    Pi1 = 1.D0

    Qext = 0.D0
    Qsca = 0.D0
    dQextdx1 = 0.D0
    dQextdx2 = 0.D0
    dQextdRem = 0.D0
    dQextdImm = 0.D0
    dQscadx1 = 0.D0
    dQscadx2 = 0.D0
    dQscadRem = 0.D0
    dQscadImm = 0.D0

    if ARG_PRESENT(asym) or $
        ARG_PRESENT(dasymdx) or $
        ARG_PRESENT(dasymdRem) or $
        ARG_PRESENT(dasymdImm) then begin

        asym = 0.D0
        dgdx = 0.D0
        dgdn = 0.D0
        dgdk = 0.d0

        calc_g = 1b
    endif else calc_g = 0b

    Sp = 0.D0
    Sm = 0.D0
    dSpdx = 0.D0
    dSmdx = 0.D0
    dSpdRem = 0.D0
    dSpdImm = 0.D0
    dSmdRem = 0.D0
    dSmdImm = 0.D0

    for n = 1,Nstop do begin
        psi1 = double(2d0*n-1d0)*psi0/x - psim1
        chi1 = double(2d0*n-1d0)*chi0/x - chim1 ;the recurrence relations

        zeta = dcomplex(psi1,chi1)
        zetanm1 = dcomplex(psi0,chi0)

        a_n = ((D(n)/m+n/x)*psi1-psi0) / ((D(n)/m+n/x)*zeta-zetanm1)
        b_n = ((D(n)*m+n/x)*psi1-psi0) / ((D(n)*m+n/x)*zeta-zetanm1) ; formulae D9,D10

        i = dcomplex(0,1)
        dadx = i*((D(n)^2d0)*(1/m^2d0-1d0) + n*(n+1d0)*(1d0/y^2d0-1d0/x^2d0)) / $
               ((D(n)/m+n/x)*zeta-zetanm1)^2d0
        dbdx = i*(1- m^2d0 + n*(n+1d0)*(m^2d0/y^2d0-1d0/x^2d0))               / $
               ((D(n)*m+n/x)*zeta-zetanm1)^2d0

        dadRem = i*(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) - D(n)) / $
                 ((D(n)+n*m/x)*zeta - m*zetanm1)^2d0
        dadImm =  -(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) - D(n)) / $
                 ((D(n)+n*m/x)*zeta - m*zetanm1)^2d0
        dbdRem = i*(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) + D(n)) / $
                 ((D(n)*m+n/x)*zeta-zetanm1)^2d0
        dbdImm =  -(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) + D(n)) / $
                 ((D(n)*m+n/x)*zeta-zetanm1)^2d0 ; my formulae for the derivatives

        Qext = (2d0*n+1d0)*double(a_n + b_n) + Qext

        dQextdx1  = (2d0*n+1d0)*double(dadx + dbdx)    + dQextdx1
        dQextdx2  = (2d0*n+1d0)*double(a_n + b_n)      + dQextdx2
        dQextdRem = (2d0*n+1d0)*double(dadRem +dbdRem) + dQextdRem
        dQextdImm = (2d0*n+1d0)*double(dadImm +dbdImm) + dQextdImm

        Qsca = (2*n+1)*double(a_n*conj(a_n) + b_n*conj(b_n)) + Qsca

        dQscadx1 = (2d0*n+1d0)*(double(a_n)*double(dadx) + $
                   imaginary(a_n)*imaginary(dadx) + double(b_n)*double(dbdx) + $
                   imaginary(b_n)*imaginary(dbdx)) + dQscadx1

        dQscadx2 = (2d0*n+1d0)*double(a_n*conj(a_n) + b_n*conj(b_n)) + dQscadx2

        dQscadRem = (2d0*n+1d0)*(double(a_n)*double(dadRem) + $
                    imaginary(a_n)*imaginary(dadRem) + double(b_n)*double(dbdRem) + $
                    imaginary(b_n)*imaginary(dbdRem)) + dQscadRem

        dQscadImm = (2d0*n+1d0)*(double(a_n)*double(dadImm) + $
                    imaginary(a_n)*imaginary(dadImm) + double(b_n)*double(dbdImm) + $
                    imaginary(b_n)*imaginary(dbdImm)) + dQscadImm

        if calc_g then begin
           if (n gt 1) then begin
              ; We are running one n behind for asymmetry parameter, since
              ; we require a_(n+1), b_(n+1) etc...
              dn = double( n-1 )

              pre0 = dn * ( 2d0 + dn ) / (dn + 1d0)
              pre1 = ( 2d0*dn + 1d0) / dn / (dn + 1d0)


              ; At the end, we'll do asym=asym * 4/x^2
              asym += pre0 * double( anm1*conj(a_n) + bnm1*conj(b_n) ) + $
                      pre1 * double( anm1 * conj(bnm1) )


              ; Clean this up at the end. These are not the final expressions!
              dgdx += pre0 * double( conj(a_n)*dadxm1 + anm1*conj(dadx) +$
                                     conj(b_n)*dbdxm1 + bnm1*conj(dbdx)   ) + $
                      pre1 * double( dadxm1*conj(bnm1) + anm1*conj(dbdxm1) )

              dgdn += pre0 * double( conj(a_n)*dadnm1 + anm1*conj(dadRem) +$
                                     conj(b_n)*dbdnm1 + bnm1*conj(dbdRem)   ) + $
                      pre1 * double( dadnm1*conj(bnm1) + anm1*conj(dbdnm1) )

              dgdk += pre0 * double( conj(a_n)*dadkm1 + anm1*conj(dadImm) +$
                                     conj(b_n)*dbdkm1 + bnm1*conj(dbdImm)   ) + $
                      pre1 * double( dadkm1*conj(bnm1) + anm1*conj(dbdkm1) )
           endif

           anm1 = a_n
           bnm1 = b_n
           dadxm1 = dadx
           dbdxm1 = dbdx
           dadnm1 = dadRem
           dbdnm1 = dbdRem
           dadkm1 = dadImm
           dbdkm1 = dbdImm
        endif

        ; formulae for the derivatives of the Qs

        psim1 = psi0
        psi0 = psi1
        chim1 = chi0
        chi0 = chi1

        s = Dqv*Pi1
        t = s - Pi0
        Tau1 = n*t - Pi0
        Sp      = (2d0*n+1d0)*(a_n+b_n)*(Pi1+Tau1)/(n^2d0+n) + Sp
        Sm      = (2d0*n+1d0)*(a_n-b_n)*(Pi1-Tau1)/(n^2d0+n) + Sm
        dSpdx   = (2d0*n+1d0)*(dadx+dbdx)*(Pi1+Tau1)/(n^2d0+n) + dSpdx
        dSmdx   = (2d0*n+1d0)*(dadx-dbdx)*(Pi1-Tau1)/(n^2d0+n) + dSmdx
        dSpdRem = (2d0*n+1d0)*(dadRem+dbdRem)*(Pi1+Tau1)/(n^2d0+n) + dSpdRem
        dSpdImm = (2d0*n+1d0)*(dadImm+dbdImm)*(Pi1+Tau1)/(n^2d0+n) + dSpdImm
        dSmdRem = (2d0*n+1d0)*(dadRem-dbdRem)*(Pi1-Tau1)/(n^2d0+n) + dSmdRem
        dSmdImm = (2d0*n+1d0)*(dadImm-dbdImm)*(Pi1-Tau1)/(n^2d0+n) + dSmdImm
        Pi0 = Pi1
        Pi1 = s + t*(n+1d0)/n
    endfor

    Qext = 2d0*Qext/x^2d0

    dQextdx1 = 2d0*dQextdx1/x^2d0
    dQextdx2 = 4d0*dQextdx2/x^3d0
    dQextdx = dQextdx1 - dQextdx2

    dQextdRem = 2d0*dQextdRem/x^2d0
    dQextdImm = 2d0*dQextdImm/x^2d0

    Qsca = 2d0*Qsca/x^2d0

    dQscadx1 = 4d0*dQscadx1/x^2d0
    dQscadx2 = 4d0*dQscadx2/x^3d0
    dQscadx = dQscadx1 - dQscadx2

    dQscadRem = 4d0*dQscadRem/x^2d0
    dQscadImm = 4d0*dQscadImm/x^2d0

    S1 = (Sp + Sm)/2d0
    S2 = (Sp - Sm)/2d0
    dS1dx = (dSpdx + dSmdx)/2d0
    dS2dx = (dSpdx - dSmdx)/2d0
    dS1dRem = (dSpdRem + dSmdRem)/2d0
    dS1dImm = (dSpdImm + dSmdImm)/2d0
    dS2dRem = (dSpdRem - dSmdRem)/2d0
    dS2dImm = (dSpdImm - dSmdImm)/2d0

    i1 = double(S1*conj(S1))
    i2 = double(S2*conj(S2))
    di1dx = 2*(double(S1)*double(dS1dx) + imaginary(S1)*imaginary(dS1dx))
    di2dx = 2*(double(S2)*double(dS2dx) + imaginary(S2)*imaginary(dS2dx))
    di1dRem = 2*(double(S1)*double(dS1dRem) + imaginary(S1)*imaginary(dS1dRem))
    di1dImm = 2*(double(S1)*double(dS1dImm) + imaginary(S1)*imaginary(dS1dImm))
    di2dRem = 2*(double(S2)*double(dS2dRem) + imaginary(S2)*imaginary(dS2dRem))
    di2dImm = 2*(double(S2)*double(dS2dImm) + imaginary(S2)*imaginary(dS2dImm))

    if calc_g then begin
        fxx = 4d0 / x / x

        asym = fxx * asym / Qsca

        dasymdx   = ( fxx*dgdx - asym*(dQscadx+2d0*Qsca/x) ) / Qsca
        dasymdRem = ( fxx*dgdn - asym*dQscadRem ) / Qsca
        dasymdImm = ( fxx*dgdk - asym*dQscadImm ) / Qsca
    endif

    return
end
