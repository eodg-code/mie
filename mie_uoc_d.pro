pro mie_uoc_d, Dx,Cm,Inp,Dqv,Xs1,Xs2,Dqxt,Dqsc,Dqbk,Dg,Dph
;+
; NAME:
;     MIE_UOC_D
;
; PURPOSE:
;     Calculates the scattering parameters of a series of particles
;     using the Mie scattering theory.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;     mie_uoc_d, Dx, Cm, Inp, Dqv $
;     [, Xs1] [, Xs2] [, Dqxt] [, Dqsc] [, Dqbk] [, Dg] [, Dph]
;
; INPUTS:
;     Dx:   A 1D array of particle size parameters
;     Cm:   The complex refractive index of the partilces
;     Inp:  Number of scattering angles at which to calculate intensity
;           functions etc
;     Dqv:  The cosine of the scattering angles at which to calculate
;           the intensity functions etc
;
; KEYWORD INPUTS:
;
; OUTPUTS:
;     Xs1:  The first amplitude function - amplitude of light polarised
;           in the plane perpendicular to the directions of inicident
;           light propogation and observation.
;     Xs2:  The second amplitude function - amplitude of light polarised
;           in the plane parallel to the directions of inicident light
;           propogation and observation. NB. Xs1 and Xs2 are complex
;           arrays of the same dimension as Dqv
;     Dqxt: The extinction efficiency
;     Dqsc: The scattering efficiency
;     Dg:   The asymetry parameter
;     Dph:  The phase function - an array of the same dimension as Dqv.
;
; KEYWORD OUTPUTS:
;
; MODIFICATION HISTORY:
;       G. Thomas, 1998: mie_uoc.pro (translation of mieint.f to IDL)
;       D. Grainger, 2001: mie_uoc_d.pro (Added support for arrays of
;           particle sizes and included calculation of phase function)
;       G. Thomas, Sept 2003: (Put into EODG routines format)
;-

  Imaxx = 12000
  RIMax = 2.5
  Itermax = Imaxx * RIMax
  Imaxnp = 1100 ; Change this as required
  Sizes = N_Elements(Dx)
  if (N_Elements(Inp) Eq 0) then begin
    Inp = 1
    Dqv = 0
  endif

  Xs1 = ComplexArr(Inp,Sizes)
  Xs2 = ComplexArr(Inp,Sizes)
  Dqxt = Dblarr(Sizes)
  Dqsc = Dblarr(Sizes)
  Dqbk = Dblarr(Sizes)
  Dg = Dblarr(Sizes)
  Dph = Dblarr(Inp,Sizes)

  for Size = 0, Sizes - 1 do begin

    if (Dx(Size) GT Imaxx) then message, 'Error: Size Parameter Overflow in Mie'
    Ir = 1.D0 / Cm
    Y =  Dx(Size) * Cm

    if (Dx(Size) LT 0.02) then NStop = 2 else begin
      case 1 OF
        (Dx(Size) LE 8.0)    : NStop = Dx(Size) + 4.00*Dx(Size)^(1./3.) + 2.0
        (Dx(Size) LT 4200.0) : NStop = Dx(Size) + 4.05*Dx(Size)^(1./3.) + 2.0
        else                 : NStop = Dx(Size) + 4.00*Dx(Size)^(1./3.) + 2.0
      endcase
    end
    NmX = FIX(MAX([NStop,ABS(Y)]) + 15.)
    D = DCOMPLEXARR(Nmx+1)

    for N = Nmx-1,1,-1 do begin
      A1 = (N+1) / Y
      D(N) = A1 - 1/(A1+D(N+1))
    end

    Sm = DCOMPLEXARR(Inp)
    Sp = DCOMPLEXARR(Inp)
    Pi0 = DCOMPLEXARR(Inp)
    Pi1 = DCOMPLEX(REPLICATE(1.D0,Inp),REPLICATE(0.D0,Inp))

    Psi0 = Cos(Dx(Size))
    Psi1 = Sin(Dx(Size))
    Chi0 =-Sin(Dx(Size))
    Chi1 = Cos(Dx(Size))
    Xi0 = DCOMPLEX(Psi0,Chi0)
    Xi1 = DCOMPLEX(Psi1,Chi1)

    Dg(Size) = 0.D0
    Dqsc(Size) = 0.D0
    Dqxt(Size) = 0.D0
    Tnp1 = 1

    for N = 1,Nstop do begin
      DN = Double(N)
      Tnp1 = Tnp1 + 2
      Tnm1 = Tnp1 - 2
      A2 = Tnp1 / (DN*(DN+1.D0))
      Turbo = (DN+1.D0) / DN
      Rnx = DN/Dx(Size)
      Psi = DOUBLE(Tnm1)*Psi1/Dx(Size) - Psi0
      Chi = Tnm1*Chi1/Dx(Size)       - Chi0
      Xi = DCOMPLEX(Psi,Chi)
      A = ((D[N]*Ir+Rnx)*Psi-Psi1) / ((D[N]*Ir+Rnx)*  Xi-  Xi1)
      B = ((D[N]*Cm+Rnx)*Psi-Psi1) / ((D[N]*Cm+Rnx)*  Xi-  Xi1)
      Dqxt(Size) = Tnp1 *      DOUBLE(A + B)          + Dqxt(Size)
      Dqsc(Size) = Tnp1 * DOUBLE(A*CONJ(A) + B*CONJ(B)) + Dqsc(Size)
      if (N GT 1) then Dg(Size) = Dg(Size) $
                            + (dN*dN - 1) * DOUBLE(ANM1*CONJ(A) + BNM1 * CONJ(B)) / dN $
                            + TNM1 * DOUBLE(ANM1*CONJ(BNM1)) / (dN*dN - dN)
      Anm1 = A
      Bnm1 = B

      S = Dqv * Pi1
      T = S - Pi0
      Taun = N*T - Pi0
      Sp = (A2 * (A + B)) * (Pi1 + Taun) + Sp
      Sm = (A2 * (A - B)) * (Pi1 - Taun) + Sm
      Pi0 = Pi1
      Pi1 = S + T*Turbo

      Psi0 = Psi1
      Psi1 = Psi
      Chi0 = Chi1
      Chi1 = Chi
      Xi1 = DCOMPLEX(Psi1,Chi1)

    end ; for Nstop

    if (Dg(Size) GT 0) then Dg(Size) = 2 * Dg(Size) / Dqsc(Size)
    Dqsc(Size) =  2 * Dqsc(Size) / Dx(Size)^2
    Dqxt(Size) =  2 * Dqxt(Size) / Dx(Size)^2

    Xs1(*,Size) = (Sp + Sm) / 2
    Xs2(*,Size) = (Sp - Sm) / 2
    Dph(*,Size) = 2 * DOUBLE(Xs1(*,Size)*CONJ(Xs1(*,Size)) + Xs2(*,Size)*CONJ(Xs2(*,Size))) / (Dx(Size)^2 * Dqsc(Size))
    Dqbk(Size) =  4 * ABS(Xs1(Inp-1)^2) / Dx(Size)^2
  end

  return
end
