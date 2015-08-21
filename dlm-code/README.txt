This directory contains the DLM version of the mie_single routine.
The object library and dlm file for both x86 (32 bit) and x86_64 (64 
bit) linux are mirrored in the directory. See the instructions in that
directory to use the code.
/home/crun/eodg/idl/dlm/

The DLM routine can be called in the same way as the normal 
mie_single routine:
mie_dlm_single, Dx, Cm, [Dqv=dqv], Dqxt, Dqsc, Dqbk, Dg, Xs1, Xs2, Dph

Dx must be a vector of type double (it can have one element though)
Cm must be a scaler of type double complex
Dqv (if it is specified) must be a vector of type double (it can have 
one element though)

Note that this has only been compiled for use on Linux, but it should be 
possible to compile the code for use with the Alpha systems without too 
many alterations. In fact it should be even work on Windows machines 
(with a few alterations), you are so perverse as to want to use IDL 
under Windows.
