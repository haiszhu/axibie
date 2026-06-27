# AxiStokes3D build Makefile.  Compiles the Fortran in src/ + the mwrap MEX gateway and links
# the support libraries build/libaxissymstok_kernelsplit.a and build/libaxissymlap_kernelsplit.a,
# which are shipped prebuilt together with their module interfaces
# build/axissymstok_kernelsplit_mod.mod and build/axissymlap_kernelsplit_mod.mod (the .mod lets
# the sources that `use` those modules compile; the .a provides the compiled routines at link).

PROJECT := AxiStokes3D
PP      := ax

FC := gfortran-15
CC := gcc-15
AR := ar
MW := $(HOME)/mwrap/mwrap
MWFLAGS := -c99complex -i8 -mex

SRC_DIR    := src
MATLAB_DIR := matlab
TEST_DIR   := test
BLD_DIR    := build

FFLAGS  := -O3 -fPIC -fdefault-integer-8 -fallow-argument-mismatch \
           -frecursive -cpp -ffp-contract=off -fno-unsafe-math-optimizations \
           -march=native -funroll-loops -fopenmp \
           -J$(BLD_DIR) -I$(BLD_DIR) -I$(SRC_DIR)

UNAME := $(shell uname)
ARCH  := $(shell uname -m)
ifeq ($(UNAME), Darwin)
  MATLAB_ROOT := $(shell ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort | tail -n1)
  ifeq ($(ARCH), arm64)
    MATLAB_ARCH  := maca64
    MEX_EXT      := mexmaca64
    OPENBLAS_DIR := /opt/homebrew/opt/openblas-singlethread
  else
    MATLAB_ARCH  := maci64
    MEX_EXT      := mexmaci64
    OPENBLAS_DIR := /usr/local/opt/openblas-singlethread
  endif
  MATLAB_LIBS := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmx.dylib \
                 $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmex.dylib \
                 $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmat.dylib -lm
  LAPACK_LIBS := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmwlapack.dylib $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmwblas.dylib
  MEX_LDFLAGS := -bundle -Wl,-undefined,dynamic_lookup
else
  MATLAB_ROOT := $(shell ls -d /usr/local/MATLAB/R* 2>/dev/null | sort | tail -n1)
  MATLAB_ARCH := glnxa64
  MEX_EXT     := mexa64
  OPENBLAS_DIR := /usr/local/opt/openblas-singlethread
  MATLAB_LIBS := -L$(MATLAB_ROOT)/bin/$(MATLAB_ARCH) -lmx -lmex -lmat -lm
  LAPACK_LIBS := -L$(MATLAB_ROOT)/bin/$(MATLAB_ARCH) -lmwlapack -lmwblas
  MEX_LDFLAGS := -shared
endif

OPENBLAS_LIBS := -L$(OPENBLAS_DIR)/lib -lopenblas

# Fortran sources compiled here.  The axissymstok / axissymlap specialquad + kernelsplit_mex
# wrappers use their kernelsplit module -> they compile against the shipped
# build/axissym{stok,lap}_kernelsplit_mod.mod and link build/libaxissym{stok,lap}_kernelsplit.a.
PUBLIC_SOURCES := $(SRC_DIR)/axistokes3d_mod.f90 \
                  $(SRC_DIR)/axisym_modal_green_mod.f90 \
                  $(SRC_DIR)/axilaplace3d_mod.f90 \
                  $(SRC_DIR)/specialquad_mod.f90 \
                  $(SRC_DIR)/axissymstok_specialquad_mod.f90 \
                  $(SRC_DIR)/axissymstok_specialquad_mex.f90 \
                  $(SRC_DIR)/axissymstok_kernelsplit_mex.f90 \
                  $(SRC_DIR)/axissymlap_kernelsplit_mex.f90 \
                  $(SRC_DIR)/axissymlap_specialquad_mod.f90 \
                  $(SRC_DIR)/axissymlap_specialquad_mex.f90 \
                  $(SRC_DIR)/axissym_physop_mod.f90 \
                  $(SRC_DIR)/axissym_physop_mex.f90
PUBLIC_OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BLD_DIR)/%.o,$(PUBLIC_SOURCES))

PUB_LIB     := $(BLD_DIR)/lib$(PROJECT).a
SEC_LIB     := $(BLD_DIR)/libaxissymstok_kernelsplit.a   # shipped prebuilt in build/ — not a build target
SEC_LIB_LAP := $(BLD_DIR)/libaxissymlap_kernelsplit.a    # shipped prebuilt in build/ — not a build target

MW_SRC  := $(MATLAB_DIR)/$(PROJECT).mw
MEX_C   := $(MATLAB_DIR)/$(PROJECT)_mex.c
MEX_OUT := $(MATLAB_DIR)/$(PROJECT)_mex.$(MEX_EXT)

.PHONY: all lib mex clean

all:
	@echo "$(PROJECT) build targets: make mex | make lib | make clean"
	@echo "  NOTE: build/libaxissym{stok,lap}_kernelsplit.a + axissym{stok,lap}_kernelsplit_mod.mod must be present."

mex: $(MEX_OUT)
lib: $(PUB_LIB)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

$(BLD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BLD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Build-order deps (axissymstok_kernelsplit_mod.mod is shipped in $(BLD_DIR), so the files that
# `use` that module need no .o dependency -- it is already present).
$(BLD_DIR)/axisym_modal_green_mod.o: $(BLD_DIR)/axistokes3d_mod.o
$(BLD_DIR)/specialquad_mod.o: $(BLD_DIR)/axistokes3d_mod.o
$(BLD_DIR)/axissymstok_specialquad_mod.o: $(BLD_DIR)/axistokes3d_mod.o $(BLD_DIR)/specialquad_mod.o $(BLD_DIR)/axisym_modal_green_mod.o
$(BLD_DIR)/axissymstok_specialquad_mex.o: $(BLD_DIR)/axissymstok_specialquad_mod.o
$(BLD_DIR)/axilaplace3d_mod.o: $(BLD_DIR)/axistokes3d_mod.o
$(BLD_DIR)/axissymlap_specialquad_mod.o: $(BLD_DIR)/axilaplace3d_mod.o $(BLD_DIR)/specialquad_mod.o $(BLD_DIR)/axisym_modal_green_mod.o
$(BLD_DIR)/axissymlap_specialquad_mex.o: $(BLD_DIR)/axissymlap_specialquad_mod.o
$(BLD_DIR)/axissym_physop_mod.o: $(BLD_DIR)/axistokes3d_mod.o $(BLD_DIR)/axissymlap_specialquad_mod.o $(BLD_DIR)/axissymstok_specialquad_mod.o
$(BLD_DIR)/axissym_physop_mex.o: $(BLD_DIR)/axissym_physop_mod.o

$(PUB_LIB): $(PUBLIC_OBJECTS)
	$(AR) rcs $@ $^

$(MEX_C): $(MW_SRC) | $(BLD_DIR)
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) $(PROJECT)_mex -mb -list $(PROJECT).mw
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) $(PROJECT)_mex -c $(PROJECT)_mex.c $(PROJECT).mw
	perl -pi -e 's/_{2,}/_/g' $(MEX_C)

# Links the prebuilt support lib (must exist in $(BLD_DIR)); make errors loudly if it is missing.
$(MEX_OUT): $(PUB_LIB) $(SEC_LIB) $(SEC_LIB_LAP) $(MEX_C)
	$(CC) $(MEX_LDFLAGS) -fPIC \
	  -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2018a -DMX_COMPAT_32=0 \
	  -I$(MATLAB_ROOT)/extern/include \
	  $(MEX_C) \
	  -L$(BLD_DIR) -l$(PROJECT) -laxissymstok_kernelsplit -laxissymlap_kernelsplit \
	  $(MATLAB_LIBS) $(LAPACK_LIBS) \
	  -lgfortran -lm \
	  -o $(MEX_OUT)

clean:
	rm -rf $(BLD_DIR) $(MEX_OUT)
