OVERVIEW
--------

This file describes the Mie routines in this directory that are written in IDL
or FORTRAN (but called from IDL).

The routines available are:

Mie Scattering routines:
    See the mie routines web page for more information
    (www.atm.ox.ac.uk/code/mie).
mie_single.pro
    Calculates scattering parameters for a single particle.
mie_derivs.pro
    Calculates scattering parameters and derivatives for a single particle.
mie_lognormal.pro
    Calculates scattering parameters for a modified gamma or log-normal
    distribution of particles.
mie_derivs_ln.pro
    Calculates scattering parameters and derivatives for a modified gamma or
    log-normal distribution of particles.
DLM/mie_single.dlm etc
    A IDL DLM version of mie_single.pro. NOTE: This procedure no longer returns
    the phase function, as this was causing segmentation faults. The scattering
    functions Xs1 and Xs2 can be used to calculate the phase function outside
    the routine though.

Old Mie scattering code:

mieext.pro
    Calculates the extinction and scattering coefficients, plus the asymmetry
    parameter of a single particle.
mieext_f.pro
    Uses CALL_EXTERNAL to call precompiled versions of mieext.

mie_uoc_d.pro
    Similar to mie_single.pro. Retained for compatibility.

Legendre expansion of a (phase) function:
legcrt.pro
    Reconstructs a function given a set of Legendre coefficients.
legexp.pro
    Expands a function as a Legendre series, returning the Legendre coefficients
quadrature.pro
    Assigns abscissa and weights for integration on the interval [-1,1].
    Available quadrature types are:
        Simpson's
        Trapezium
        Gaussian
        Radau
        Lobarto
Shift_quadrature.pro
    Shifts the quadrature interval from [-1,1] to an arbitrary [a,b].


COMPILATION
-----------
cp make.inc.example make.inc

Edit make.inc appropriately.

make

Note: The default in make.inc.example is for the EODG computing environment.
