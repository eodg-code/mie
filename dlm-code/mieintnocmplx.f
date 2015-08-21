!     General purpose Mie scattering routine for single particles
!     Author: R Grainger 1990
!     History:
!     G Thomas, March 2005: Added calculation of Phase function and
!     code to ensure correct calculation of backscatter coeficient
!     G Thomas, November 2006: Calculation of backscatter efficiency 
!     now done using the A and B values rather than the intensity 
!     at 180 degrees.

      Options/Extend_Source
      Subroutine MieIntNoCmplx(Dx, SCm, Inp, Dqv, Dqxt, Dqsc, Dbsc, Dg, Xs1r, Xs1i, Xs2r, Xs2i, DPh, Error)

      Integer * 2  Imaxx
      Parameter (Imaxx = 12000)
      Real * 4     RIMax          ! largest real part of refractive index
      Parameter (RIMax = 2.5)
      Real * 4     IRIMax         ! largest imaginary part of refractive index
      Parameter (IRIMax = -2)
      Integer * 2  Itermax
      Parameter (Itermax = 12000 * 2.5)
                                ! must be large enough to cope with the
                                ! largest possible nmx = x * abs(scm) + 15
                                ! or nmx =  Dx + 4.05*Dx**(1./3.) + 2.0
      Integer * 2  Imaxnp
      Parameter (Imaxnp = 10000)  ! Change this as required
*     INPUT
      Real * 8     Dx
      Complex * 16  SCm
      Integer * 4  Inp
      Real * 8     Dqv(Inp)
*     OUTPUT
      Real * 8     Xs1r(Inp)
      Real * 8     Xs1i(Inp)
      Real * 8     Xs2r(Inp)
      Real * 8     Xs2i(Inp)
      Real * 8     Dqxt
      Real * 8     Dqsc
      Real * 8     Dg
      Real * 8     Dbsc
      Real * 8     DPh(InP)
      Integer * 4  Error
*     LOCAL
      Integer * 2  I
      Integer * 2  NStop
      Integer * 2  NmX
      Integer * 4  N    ! N*N > 32767 ie N > 181
      Integer * 4  Inp2
      Real * 8     Chi,Chi0,Chi1
      Real * 8     APsi,APsi0,APsi1
      Real * 8     Pi0(Imaxnp)
      Real * 8     Pi1(Imaxnp)
      Real * 8     Taun(Imaxnp)
      Real * 16    Psi,Psi0,Psi1
      Complex * 8  Ir
      Complex * 16 Xs1(InP)
      Complex * 16 Xs2(InP)
      Complex * 16 Cm
      Complex * 16 A,ANM1,APB
      Complex * 16 B,BNM1,AMB
      Complex * 16 D(Itermax)
      Complex * 16 Sp(Imaxnp)
      Complex * 16 Sm(Imaxnp)
      Complex * 16 Xi,Xi0,Xi1
      Complex * 16 Bscnum
      Complex * 16 Y
*     ACCELERATOR VARIABLES
      Integer * 2  Tnp1
      Integer * 2  Tnm1
      Real * 16    Dn
      Real * 8     Rnx
      Real * 8     S(Imaxnp)
      Real * 8     T(Imaxnp)
      Real * 8     Turbo
      Real * 8     A2
      Complex * 16 A1

      If ((Dx.Gt.Imaxx) .Or. (InP.Gt.ImaxNP)) Then
        Error = 1
        Return
      EndIf
      Cm = SCm
      Ir = 1 / Cm
      Y =  Dx * Cm
      If (Dx.Lt.0.02) Then
         NStop = 2
      Else
         If (Dx.Le.8.0) Then
            NStop = Dx + 4.00*Dx**(1./3.) + 2.0
         Else
            If (Dx.Lt. 4200.0) Then
               NStop = Dx + 4.05*Dx**(1./3.) + 2.0
            Else
               NStop = Dx + 4.00*Dx**(1./3.) + 2.0
            End If
         End If
      End If
      NmX = Max(Real(NStop),Real(Abs(Y))) + 15.
      If (Nmx .gt. Itermax) then
          Error = 1
          Return
      End If
      Inp2 = Inp+1
      D(NmX) = Dcmplx(0,0)
      Do N = Nmx-1,1,-1
         A1 = (N+1) / Y
         D(N) = A1 - 1/(A1+D(N+1))
      End Do
      Do I =1,Inp2
         Sm(I) = Dcmplx(0,0)
         Sp(I) = Dcmplx(0,0)
         Pi0(I) = 0
         Pi1(I) = 1
      End Do
      Psi0 = Cos(Dx)
      Psi1 = Sin(Dx)
      Chi0 =-Sin(Dx)
      Chi1 = Cos(Dx)
      APsi0 = Psi0
      APsi1 = Psi1
      Xi0 = Dcmplx(APsi0,Chi0)
      Xi1 = Dcmplx(APsi1,Chi1)
      Dg = 0
      Dqsc = 0
      Dqxt = 0
      Bscnum = Dcmplx(0,0)
      Tnp1 = 1
      Do N = 1,Nstop
         DN = N
         Tnp1 = Tnp1 + 2
         Tnm1 = Tnp1 - 2
         A2 = Tnp1 / (DN*(DN+1D0))
         Turbo = (DN+1D0) / DN
         Rnx = DN/Dx
         Psi = Dble(Tnm1)*Psi1/Dx - Psi0
         APsi = Psi
         Chi = Tnm1*Chi1/Dx       - Chi0
         Xi = Dcmplx(APsi,Chi)
         A = ((D(N)*Ir+Rnx)*APsi-APsi1) / ((D(N)*Ir+Rnx)*  Xi-  Xi1)
         B = ((D(N)*Cm+Rnx)*APsi-APsi1) / ((D(N)*Cm+Rnx)*  Xi-  Xi1)
         Dqxt   = Tnp1 *      Dble(A + B)          + Dqxt
         Dqsc   = Tnp1 * (A*Conjg(A) + B*Conjg(B)) + Dqsc
         Bscnum = Tnp1 * (-1**N) * (B - A)         + Bscnum
         If (N.Gt.1) then
            Dg = Dg + (dN*dN - 1) * Dble(ANM1*Conjg(A) + BNM1 * Conjg(B)) / dN + TNM1 * Dble(ANM1*Conjg(BNM1)) / (dN*dN - dN)
         End If
         Anm1 = A
         Bnm1 = B
         APB = A2 * (A + B)
         AMB = A2 * (A - B)
         Do I = 1,Inp2
            If (I.GT.Inp) Then
               S(I) = -Pi1(I)
            Else
               S(I) = Dqv(I) * Pi1(I)
            End If
            T(I) = S(I) - Pi0(I)
            Taun(I) = N*T(I) - Pi0(I)
            Sp(I) = APB * (Pi1(I) + Taun(I)) + Sp(I)
            Sm(I) = AMB * (Pi1(I) - Taun(I)) + Sm(I)
            Pi0(I) = Pi1(I)
            Pi1(I) = S(I) + T(I)*Turbo
         End Do
         Psi0 = Psi1
         Psi1 = Psi
         Apsi1 = Psi1
         Chi0 = Chi1
         Chi1 = Chi
         Xi1 = Dcmplx(APsi1,Chi1)
      End Do
      If (Dg .GT.0) Dg = 2 * Dg / Dqsc
      Dqsc =  2 * Dqsc / Dx**2
      Dqxt =  2 * Dqxt / Dx**2
      Dbsc =  Dble(Bscnum*Conjg(Bscnum)) / Dx**2
      Do I = 1,Inp
         Xs1(I) = (Sp(I)+Sm(I)) / 2
         Xs2(I) = (Sp(I)-Sm(I)) / 2
         Dph(I) = 2 * Dble(Xs1(I)*Conjg(Xs1(I)) + Xs2(I)*Conjg(Xs2(I))) / (Dx**2 * Dqsc)
         Xs1r(I) = Real(Xs1(I))
         Xs1i(I) = AImag(Xs1(I))
         Xs2r(I) = Real(Xs2(I))
         Xs2i(I) = AImag(Xs2(I))
      End Do
      Error = 0
      Return
      End
