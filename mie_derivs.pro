PRO mie_derivs, x,Cm,Dqv,Qext,Qsca,$
                dQextdx,dQextdRem,dQextdImm,$
                dQscadx,dQscadRem,dQscadImm,$
                i1,i2,di1dx,di2dx,di1dRem,di1dImm,di2dRem,di2dImm,$
                silent=silent, $
                asym=asym,dasymdx=dasymdx, $
                dasymdRem=dasymdRem,dasymdImm=dasymdImm

;+
; NAME:
;     MIE_DERIVS
;
; PURPOSE:
;     Calculates the scattering parameters and their analytical
;     derivatives of a series of particles using Mie scattering theory.
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
;     mie_derivs, x, Cm, Dqv $
;     [, Qext][, Qsca][, dQextdx][, dQextdRem][, dQextdImm] $
;     [, dQscadx][, dQscadRem][, dQscadImm][, i1][, i2] $
;     [, di1dx][, di2dx][, di1dRem][, di1dImm][, di2dRem][, di2dImm] $
;     [, ASYM=asym][, dASYMdx=dasymdx][, dASYMdRem=dasymdrem] $
;     [, dASYMdImm=dASYMdrem][, /SILENT]
;
; INPUTS:
;     x:         The particle size parameter
;     Cm:        The complex refractive index of the particle
;     Dqv:       The cosine of the scattering angles at which to
;                calculate the intensity functions etc
;
; KEYWORD INPUTS:
;     SILENT:    Don't give warnings about positive k values.
;
; OUTPUTS:
;     Qext:      The extinction efficiency
;     Qsca:      The extinction efficiency
;     dQextdx:   Derivative of the extinction efficiency wrt the
;                particle size parameter
;     dQextdRem: Derivative of the extinction efficiency wrt the
;                real part of the refractive index
;     dQextdImm: Derivative of the extinction efficiency wrt the
;                imaginary part of the refractive index
;     dQscadx:   Derivative of the scattering efficiency wrt the
;                particle size parameter
;     dQscadRem: Derivative of the scattering efficiency wrt the
;                real part of the refractive index
;     dQscadImm: Derivative of the scattering efficiency wrt the
;                imaginary part of the refractive index
;     i1:        The first intensity function - intensity of light
;                polarized in the plane perpendicular to the directions
;                of incident light propagation and observation.
;     i2:        The second intensity function - intensity of light
;                polarized in the plane parallel to the directions of
;                incident light propagation and observation.
;     di1dx:     Derivatives of the intensity functions wrt the
;     di2dx:     particle size parameter
;     di1dRem:   Derivatives of the first intensity function wrt to
;     di1dImm:   the real and imaginary parts of the refractive index
;     di2dRem:   Derivatives of the second intensity function wrt
;     di2dImm:   to the real and imaginary parts of the refractive index
;                NB. The values of involving the intensity functions
;                are arrays of the same dimension as dqv
;
; KEYWORD OUTPUTS:
;     asym:       The asymmetry parameter.
;     dasymdx:    Derivative of asymmetry wrt size parameter.
;     dasymdRem:  Derivative of asymmetry wrt real pt of RI.
;     dasymdImm:  Derivative of asymmetry wrt imag pt of RI.
;
; RESTRICTIONS:
;
; MODIFICATION HISTORY:
;     J. Lucas, 2002: miederivs.pro (Main line programme)
;     D. Grainger, 2002: pro_miederivs.pro (Procedure version of
;         miederivs.pro)
;     G. Thomas, Sep 2003: mie_derivs.pro (Put into EODG routines
;         format)
;     J. Graham (UC Berkeley), Feb 2004: (Introduced explicit double
;         precision numerical values into all computational expressions)
;     G. Thomas, Feb 2004: Header updated.
;     A. Smith, Aug 2010: Added positive k warning.
;     A. Smith, Apr 2013: Added asymmetry parameter.
;-

if imaginary(cm) gt 0d0 and not(keyword_set(silent)) then $
    message, /continue,'Warning: Imaginary part of refractive index '+$
             'should be negative for absorbing particles. '+$
             'Set /silent to hide this message.'

IF x LT 0.02 THEN NStop = 2 ELSE BEGIN
    CASE 1 OF
        x LE 8.0    : NStop = x + 4.00*x^(1./3.) + 1.0
        x LT 4200.0 : NStop = x + 4.05*x^(1./3.) + 2.0
        ELSE        : NStop = x + 4.00*x^(1./3.) + 2.0 ;giving no. of terms required for convergence of Mie expressions (Wiscombe 1980)
    ENDCASE
ENDELSE

m = cm

y = m*x
Nmx = FIX(MAX([NStop,ABS(y)]) + 15.)
D = DCOMPLEXARR(Nmx+1)

FOR p = Nmx-1,1,-1 DO BEGIN
    R = (p+1) / y
    D(p) = R - 1/(R+D(p+1)) ;downward recurrence to find A_n(y)
ENDFOR

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

IF ARG_PRESENT(asym) OR $
   ARG_PRESENT(dasymdx) OR $
   ARG_PRESENT(dasymdRem) OR $
   ARG_PRESENT(dasymdImm) THEN BEGIN

   asym = 0.D0
   dgdx = 0.D0
   dgdn = 0.D0
   dgdk = 0.d0

   calc_g = 1b
ENDIF ELSE calc_g = 0b


Sp = 0.D0
Sm = 0.D0
dSpdx = 0.D0
dSmdx = 0.D0
dSpdRem = 0.D0
dSpdImm = 0.D0
dSmdRem = 0.D0
dSmdImm = 0.D0

FOR n = 1,Nstop DO BEGIN

    psi1 = Double(2d0*n-1d0)*psi0/x - psim1
    chi1 = Double(2d0*n-1d0)*chi0/x - chim1 ;the recurrence relations

    zeta = DCOMPLEX(psi1,chi1)
    zetanm1 = DCOMPLEX(psi0,chi0)

    a_n = ((D(n)/m+n/x)*psi1-psi0) / ((D(n)/m+n/x)*zeta-zetanm1)
    b_n = ((D(n)*m+n/x)*psi1-psi0) / ((D(n)*m+n/x)*zeta-zetanm1) ;formulae D9,D10

    i = dcomplex(0,1)
    dadx = i*((D(n)^2d0)*(1/m^2d0-1d0) + n*(n+1d0)*(1d0/y^2d0-1d0/x^2d0)) / ((D(n)/m+n/x)*zeta-zetanm1)^2d0
    dbdx = i*(1- m^2d0 + n*(n+1d0)*(m^2d0/y^2d0-1d0/x^2d0)) / ((D(n)*m+n/x)*zeta-zetanm1)^2d0

    dadRem = i*(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) - D(n)) / ((D(n)+n*m/x)*zeta - m*zetanm1)^2d0
    dadImm =  -(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) - D(n)) / ((D(n)+n*m/x)*zeta - m*zetanm1)^2d0
    dbdRem = i*(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) + D(n)) / ((D(n)*m+n/x)*zeta-zetanm1)^2d0
    dbdImm =  -(n*(n+1d0)/y - y*(1d0+(D(n)^2d0)) + D(n)) / ((D(n)*m+n/x)*zeta-zetanm1)^2d0   ;my formulae for the derivatives

    Qext = (2d0*n+1d0)*DOUBLE(a_n + b_n) + Qext

    dQextdx1  = (2d0*n+1d0)*DOUBLE(dadx + dbdx)    + dQextdx1
    dQextdx2  = (2d0*n+1d0)*DOUBLE(a_n + b_n)      + dQextdx2
    dQextdRem = (2d0*n+1d0)*DOUBLE(dadRem +dbdRem) + dQextdRem
    dQextdImm = (2d0*n+1d0)*DOUBLE(dadImm +dbdImm) + dQextdImm

    Qsca = (2*n+1)*DOUBLE(a_n*CONJ(a_n) + b_n*CONJ(b_n)) + Qsca

    dQscadx1 = (2d0*n+1d0)*(DOUBLE(a_n)*DOUBLE(dadx) + $
               IMAGINARY(a_n)*IMAGINARY(dadx) + DOUBLE(b_n)*DOUBLE(dbdx) + $
               IMAGINARY(b_n)*IMAGINARY(dbdx)) + dQscadx1

    dQscadx2 = (2d0*n+1d0)*DOUBLE(a_n*CONJ(a_n) + b_n*CONJ(b_n)) + dQscadx2

    dQscadRem = (2d0*n+1d0)*(DOUBLE(a_n)*DOUBLE(dadRem) + $
                IMAGINARY(a_n)*IMAGINARY(dadRem) + DOUBLE(b_n)*DOUBLE(dbdRem) + $
                IMAGINARY(b_n)*IMAGINARY(dbdRem)) + dQscadRem

    dQscadImm = (2d0*n+1d0)*(DOUBLE(a_n)*DOUBLE(dadImm) + $
                IMAGINARY(a_n)*IMAGINARY(dadImm) + DOUBLE(b_n)*DOUBLE(dbdImm) + $
                IMAGINARY(b_n)*IMAGINARY(dbdImm)) + dQscadImm


    IF calc_g THEN BEGIN
       IF (n GT 1) THEN BEGIN
          ; We are running one n behind for asymmetry parameter, since
          ; we require a_(n+1), b_(n+1) etc...
          dn = DOUBLE( n-1 )

          pre0 = dn * ( 2d0 + dn ) / (dn + 1d0)
          pre1 = ( 2d0*dn + 1d0) / dn / (dn + 1d0)


          ; At the end, we'll do asym=asym * 4/x^2
          asym += pre0 * DOUBLE( anm1*CONJ(a_n) + bnm1*CONJ(b_n) ) + $
                  pre1 * DOUBLE( anm1 * CONJ(bnm1) )


          ; Clean this up at the end. These are not the final expressions!
          dgdx += pre0 * DOUBLE( CONJ(a_n)*dadxm1 + anm1*CONJ(dadx) +$
                                 CONJ(b_n)*dbdxm1 + bnm1*CONJ(dbdx)   ) + $
                  pre1 * DOUBLE( dadxm1*CONJ(bnm1) + anm1*CONJ(dbdxm1) )

          dgdn += pre0 * DOUBLE( CONJ(a_n)*dadnm1 + anm1*CONJ(dadRem) +$
                                 CONJ(b_n)*dbdnm1 + bnm1*CONJ(dbdRem)   ) + $
                  pre1 * DOUBLE( dadnm1*CONJ(bnm1) + anm1*CONJ(dbdnm1) )

          dgdk += pre0 * DOUBLE( CONJ(a_n)*dadkm1 + anm1*CONJ(dadImm) +$
                                 CONJ(b_n)*dbdkm1 + bnm1*CONJ(dbdImm)   ) + $
                  pre1 * DOUBLE( dadkm1*CONJ(bnm1) + anm1*CONJ(dbdkm1) )


       ENDIF
       anm1 = a_n
       bnm1 = b_n
       dadxm1 = dadx
       dbdxm1 = dbdx
       dadnm1 = dadRem
       dbdnm1 = dbdRem
       dadkm1 = dadImm
       dbdkm1 = dbdImm

    ENDIF

   ;formulae for the derivatives of the Qs

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

ENDFOR

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

i1 = DOUBLE(S1*CONJ(S1))
i2 = DOUBLE(S2*CONJ(S2))
di1dx = 2*(DOUBLE(S1)*DOUBLE(dS1dx) + IMAGINARY(S1)*IMAGINARY(dS1dx))
di2dx = 2*(DOUBLE(S2)*DOUBLE(dS2dx) + IMAGINARY(S2)*IMAGINARY(dS2dx))
di1dRem = 2*(DOUBLE(S1)*DOUBLE(dS1dRem) + IMAGINARY(S1)*IMAGINARY(dS1dRem))
di1dImm = 2*(DOUBLE(S1)*DOUBLE(dS1dImm) + IMAGINARY(S1)*IMAGINARY(dS1dImm))
di2dRem = 2*(DOUBLE(S2)*DOUBLE(dS2dRem) + IMAGINARY(S2)*IMAGINARY(dS2dRem))
di2dImm = 2*(DOUBLE(S2)*DOUBLE(dS2dImm) + IMAGINARY(S2)*IMAGINARY(dS2dImm))


IF calc_g THEN BEGIN

   fxx = 4d0 / x / x

   asym = fxx * asym / Qsca

   dasymdx   = ( fxx*dgdx - asym*(dQscadx+2d0*Qsca/x) ) / Qsca
   dasymdRem = ( fxx*dgdn - asym*dQscadRem ) / Qsca
   dasymdImm = ( fxx*dgdk - asym*dQscadImm ) / Qsca

ENDIF

RETURN

END