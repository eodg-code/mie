pro mieext, NPts, Dx, Cm, Dqxt, Dqsc, Dg

    Imaxx = 2600
    Itermax = 3500
    Imaxnp = 1100 ; Change this as required

    Dqxt(*) = 0D0
    Dqsc(*) = 0D0
    Dg(*) = 0D0

    for I = 0, NPts - 1 do begin
        if (Dx(I) GT Imaxx) then message, 'Error: Size Parameter Overflow in Mie'
        Ir = 1.D0 / Cm
        Y =  Dx(I) * Cm

        if (Dx(I) LT 0.02) then NStop = 2 else $
          begin
            case 1 OF
            (Dx(I) LE 8)    : NStop = Dx(I) + 4.00*Dx(I)^(1./3.) + 2.0
            (Dx(I) LT 4200) : NStop = Dx(I) + 4.05*Dx(I)^(1./3.) + 2.0
            else            : NStop = Dx(I) + 4.00*Dx(I)^(1./3.) + 2.0
           	endcase
          end
        NmX = FIX(MAX([NStop,ABS(Y)]) + 15.)
        D = DCOMPLEXARR(Nmx+1)

        for N = Nmx-1,1,-1 do $
          begin
            A1 = (N+1) / Y
            D[N] = A1 - 1/(A1+D[N+1])
          end


        Psi0 = Cos(Dx(I))
        Psi1 = Sin(Dx(I))
        Chi0 =-Sin(Dx(I))
        Chi1 = Cos(Dx(I))
        Xi0 = DCOMPLEX(Psi0,Chi0)
        Xi1 = DCOMPLEX(Psi1,Chi1)

        Tnp1 = 1

        for N = 1,Nstop do $
          begin
            DN = Double(N)
            Tnp1 = Tnp1 + 2
            Tnm1 = Tnp1 - 2
            A2 = Tnp1 / (DN*(DN+1D0))
            Turbo = (DN+1D0) / DN
            Rnx = DN/Dx(I)
            Psi = DOUBLE(Tnm1)*Psi1/Dx(I) - Psi0
            Chi = Tnm1*Chi1/Dx(I)       - Chi0
            Xi = DCOMPLEX(Psi,Chi)
            A = ((D[N]*Ir+Rnx)*Psi-Psi1) / ((D[N]*Ir+Rnx)*  Xi-  Xi1)
            B = ((D[N]*Cm+Rnx)*Psi-Psi1) / ((D[N]*Cm+Rnx)*  Xi-  Xi1)
            Dqxt(I) = Tnp1 * DOUBLE(A + B) + Dqxt(I)
            Dqsc(I) = Tnp1 * DOUBLE(A*CONJ(A) + B*CONJ(B)) + Dqsc(I)
            if (N GT 1) then Dg(I) = Dg(I) $
                                 + (dN*dN - 1) * DOUBLE(ANM1*CONJ(A) + BNM1 * CONJ(B)) / dN $
                                 + TNM1 * DOUBLE(ANM1*CONJ(BNM1)) / (dN*dN - dN)
            Anm1 = A
            Bnm1 = B
            APB = A2 * (A + B)
            AMB = A2 * (A - B)


            Psi0 = Psi1
            Psi1 = Psi
            Chi0 = Chi1
            Chi1 = Chi
            Xi1 = DCOMPLEX(Psi1,Chi1)
          end; for Nstop

        if (Dg(I) GT 0) then Dg(I) = 2 * Dg(I) / Dqsc(I)
        Dqsc(I) =  2 * Dqsc(I) / Dx(I)^2
        Dqxt(I) =  2 * Dqxt(I) / Dx(I)^2
    endfor
end
