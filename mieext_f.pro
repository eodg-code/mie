pro mieext_f, NPts, Dx, Cm, Dqxt, Dqsc, Dg

; Multi architecture Mie routine.
; INPUT
; Npts, I4,       number of values of the size parameter (X)
; Dx,   R8(Npts), vector of size parameter values
; Cm,   C8,       complex refractive index (same for all size values)
;
; OUTPUT
; Dqxt, R8(Npts)  vector of extinction efficiency
; Dqsc, R8(Npts)  vector of scattering efficiency
; Dg,   R8(Npts)  vector of asymmetry parameter

    case !Version.arch of
        'x86'  : status = CALL_EXTERNAL('/home/crun/eodg/idl/mie/mieext_x86.so'  , $
                                        'mieextidl_', Npts, Dx, Cm, Dqxt, Dqsc, Dg)
        'alpha': status = CALL_EXTERNAL('/home/crun/eodg/idl/mie/mieext_alpha.so', $
                                        'mieextidl_', Npts, Dx, Cm, Dqxt, Dqsc, Dg)
        else: Mieext, NPts, Dx, Cm, Dqxt, Dqsc, Dg
    endcase
end
