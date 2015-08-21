This file describes the Mie routines in this directory that are
written in IDL or FORTRAN (but called from IDL).

The routines available are:

Mie Scattering routines:
   See the mie routines web page for more information (www.atm.ox.ac.uk/code/mie)
mie_single.pro
   calculates scattering parameters for a single particle
mie_derivs.pro
   calculates scattering parameters and derivatives for a single particle
mie_lognormal.pro
   calculates scattering parameters for a lognormal distribution of particles
mie_derivs_ln.pro
   calculates scattering parameters and derivatives for a lognormal distribution of particles
DLM/mie_single.dlm etc
   a IDL DLM version of mie_single.pro NOTE: that this procedure no longer returns the
   phase function, as this was causing segmentation faults. The scattering functions Xs1 and
   Xs2 can be used to calculate the phase function outside the routine though.

Old Mie scattering code:
mieext.pro
   calculates the extinction and scattering coefficients, plus the asymmetry parameter of a single particle
mieext_f.pro
   uses CALL_EXTERNAL to call precompiled versions of mieext (either mieext_x86.so or mieext_alpha.so, depending on the machine architecture)
Compilation on an alpha
f77 -c -extend_source mieext.f
ld -S -shared -o mieext_alpha.so mieext.o -lUfor -lfor -lFutil -lm -lots -lc
Compilation on an x86
ifc mieext.f -o mieext_x86.so -shared -w95 -Kpic -lm -132

mie_uoc_d.pro
   similar to mie_single.pro. Retained for compatibility.

Legendre expansion of a (phase) function:
legcrt.pro
   reconstructs a function given a set of Legendre coefficients
legexp.pro
   expands a function as a Legendre series, returning the Legendre coefficients
quadrature.pro
   assigns abscissa and weights for integration on the interval [-1,1]. Available quadrature types are:
      Simpson's
      Trapezium
      Gaussian
      Radau
      Lobarto
shift_quadrature.pro
   shifts the quadrature interval from [-1,1] to an arbitrary [a,b]

