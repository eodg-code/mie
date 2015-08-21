# mie_dlm_single IDL DLM makefile
# Gareth Thomas

# Locations of IDL and IFC may need to be changed
IDL_DIR	= /usr/local/PACK/idl80/idl/idl80
#IFC_DIR = /usr/local/PACK/intel/compiler70/ia32
#IFC_DIR = /usr/local/PACK/intel/compiler90
IFC_DIR = /usr/local/PACK/intel/compiler11.1.072

# If the version number of IDL is prior to 5.2 then set the
# IDL flag below to "-D__IDLPRE53__".
#IDL = -D__IDLPRE53__

# VARIABLES ################################################
PARALLEL = false
FCOMP = ifort
INC = -I$(IDL_DIR)/external

GCC_FLAGS = -pedantic -Wall -Wimplicit -fpic -O3 -D__G77__

ifeq ($(PARALLEL),true)
	OBJS = mie_dlm_single-openmp.o mieintnocmplx.o
#	GCC_FLAGS = -fopenmp $(GCC_FLAGS)
else
	OBJS = mie_dlm_single.o mieint.o
endif


OUT_FILE = mie_dlm_single.so

SO_EXT	= so


ifeq ($(FCOMP),g77)
#       g77 flags
	FC_FLAGS = -ff90 -fugly-complex -ffixed-line-length-none $(GCC_FLAGS)
# 	g77 compatable mieint.f
	MIESOURCE = mieint-real8.f
#	Linking flags
	LD_FLAGS = -shared -static -lm
else
# 	IFC flags
#	Version 7.0
#	FC_FLAGS = -D__G77__ -w95 -Kpic -lm -132 -O1
#	Version 9.0
	FC_FLAGS = -D__G77__ -w95 -fpic -132 -O3
#	mieint.f using Real * 16 variables
	MIESOURCE = mieint.f
# 	Linking flags for compatability with IFC
#	Version 7.0
#	LD_FLAGS = -shared -lm -lF90 -lCEPCF90 -L$(IFC_DIR)/lib
#	Version 9.0
	LD_FLAGS = -shared -lifcore -fopenmp -L$(IFC_DIR)/lib/intel64
endif


$(OUT_FILE): $(OBJS)
	gcc -o $(OUT_FILE) $(OBJS) $(LD_FLAGS)

mieint.o: mieint.f
	$(FCOMP) $(FC_FLAGS) -c -o mieint.o mieint.f
mieintnocmplx.o: mieintnocmplx.f
	$(FCOMP) $(FC_FLAGS) -c -o mieintnocmplx.o mieintnocmplx.f
mie_dlm_single.o: mie_dlm_single.c
	gcc $(GCC_FLAGS) $(INC) -c -o mie_dlm_single.o mie_dlm_single.c
mie_dlm_single-openmp.o: mie_dlm_single-openmp.c
	gcc -fopenmp $(GCC_FLAGS) $(INC) -c -o mie_dlm_single-openmp.o mie_dlm_single-openmp.c

clean:
	rm -f *.o core

