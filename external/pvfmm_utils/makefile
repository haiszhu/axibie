# Simple pvfmm mex build; defaults to Accelerate, optional OpenBLAS.
# Run `make` to list available tasks.
# g++-14 -O3 -std=c++17 -flax-vector-conversions -fPIC -fopenmp -include nla_wrap_MPI.h \
#   -I/opt/homebrew/opt/openblas/include \
#   -I/opt/homebrew/opt/open-mpi/include \
#   -I/Users/hzhu/git/pvfmm/_install/include/pvfmm \
#   -I/Users/hzhu/git/pvfmm/_install/examples/include \
#   -I/Applications/MATLAB_R2025b.app/extern/include \
#   -I/Applications/MATLAB_R2025b.app/simulink/include \
#   -DNDEBUG \
#   -dynamiclib -Wl,-undefined,dynamic_lookup -Wl,-headerpad_max_install_names \
#   -o nla_MPI.mexmaca64 \
#   nla_MPI.c nla_wrap_MPI.cpp \
#   /opt/homebrew/opt/openblas/lib/libopenblas.a \
#   /Users/hzhu/git/pvfmm/_install/lib/libpvfmm.a \
#   /opt/homebrew/opt/fftw/lib/libfftw3.a \
#   /opt/homebrew/opt/fftw/lib/libfftw3f.a \
#   -L/opt/homebrew/opt/open-mpi/lib -lmpi \
#   -L/opt/homebrew/opt/gcc@14/lib/gcc/14 -lgfortran -lquadmath \
#   -Wl,-rpath,/Applications/MATLAB_R2025b.app/bin/maca64 \
#   -Wl,-rpath,/Applications/MATLAB_R2025b.app/extern/bin/maca64 \
#   -Wl,-rpath,/opt/homebrew/opt/gcc@14/lib/gcc/14 \
#   -Wl,-rpath,/opt/homebrew/opt/open-mpi/lib \
#   /Applications/MATLAB_R2025b.app/extern/bin/maca64/libMatlabEngine.dylib \
#   /Applications/MATLAB_R2025b.app/bin/maca64/libMatlabDataArray.dylib \
#   /Applications/MATLAB_R2025b.app/bin/maca64/libmex.dylib \
#   /Applications/MATLAB_R2025b.app/bin/maca64/libmx.dylib

.DEFAULT_GOAL := all
# OS detection for reasonable defaults on Linux vs macOS.
UNAME_S := $(shell uname -s)
# Force compilers (make has a built-in default CC=cc, so use :=); auto-detect the
# latest installed homebrew GCC, matching AxiStokes3D's own Makefile/Makefile_kdtree
# convention (override with e.g. CC=gcc-14 CXX=g++-14 GFORTRAN=gfortran-14).
CC := $(shell ls /opt/homebrew/bin/gcc-[0-9]* 2>/dev/null | sort -V | tail -1)
CXX := $(shell ls /opt/homebrew/bin/g++-[0-9]* 2>/dev/null | sort -V | tail -1)
GFORTRAN := $(shell ls /opt/homebrew/bin/gfortran-* 2>/dev/null | sort -V | tail -1)
CXXFLAGS ?= -O3 -march=native -std=c++17 -flax-vector-conversions -fPIC -fopenmp
CFLAGS   ?= -O3 -fPIC -fopenmp
CPPFLAGS += -fPIC
CPPFLAGS += -DSCTL_HAVE_PVFMM -DSCTL_QUAD_T=__float128 -DSCTL_HAVE_MPI \
            -DSCTL_HAVE_BLAS -DSCTL_HAVE_LAPACK -DSCTL_HAVE_FFTW -DSCTL_HAVE_FFTWF -DSCTL_HAVE_FFTWL
CPPFLAGS += -include pvfmm_wrap.h
LDFLAGS  ?=

PVFMM_SRC_DIR   := $(abspath $(CURDIR)/../pvfmm)
PVFMM_PREFIX    := $(PVFMM_SRC_DIR)/_install
PVFMM_BUILD_DIR := $(PVFMM_SRC_DIR)/build-omp

CMAKE ?= cmake

ifneq (,$(findstring Linux,$(UNAME_S)))
  OPENBLAS_PREFIX ?= /usr
  FFTW_PREFIX     ?= /usr
  ifeq ($(OPENBLAS_PREFIX),/usr)
    OPENBLAS_LIBDIR := /usr/lib/x86_64-linux-gnu
  else
    OPENBLAS_LIBDIR := $(OPENBLAS_PREFIX)/lib
  endif  
  ifeq ($(FFTW_PREFIX),/usr)
    FFTW_LIBDIR := /usr/lib/x86_64-linux-gnu
  else
    FFTW_LIBDIR := $(FFTW_PREFIX)/lib
  endif
  CC := gcc
  CXX := g++
  GFORTRAN := gfortran
  CMAKE_C_COMPILER ?= gcc
  CMAKE_CXX_COMPILER ?= g++
  SHLIB_EXT       := so
  GFORTRAN_LIB    := libgfortran.so
  MEX_LINK_FLAGS  := -shared
  MATLAB_LIB_EXT  := so
  MATLAB_ARCH_DIR ?= glnxa64
  MEX_EXT ?= mexa64
  MEX_CC ?= gcc
  MEX_CXX ?= g++
  # Linux MATLAB ships Intel OpenMP (libiomp5) under sys/os/, NOT bin/glnxa64/libomp.
  # libiomp5 exports the GNU GOMP_* symbols, so the gateway's -fopenmp regions resolve
  # against it -> one OpenMP runtime (same strategy as make.inc.fmm3d).
  MATLAB_LIBOMP_REL := sys/os/$(MATLAB_ARCH_DIR)/libiomp5.$(MATLAB_LIB_EXT)
else
  OPENBLAS_PREFIX ?= /opt/homebrew/opt/openblas
  FFTW_PREFIX     ?= /opt/homebrew/opt/fftw
  CMAKE_C_COMPILER ?= $(CC)
  CMAKE_CXX_COMPILER ?= $(CXX)
  OPENBLAS_LIBDIR := $(OPENBLAS_PREFIX)/lib
  FFTW_LIBDIR     := $(FFTW_PREFIX)/lib
  SHLIB_EXT       := dylib
  GFORTRAN_LIB    := libgfortran.dylib
  MEX_LINK_FLAGS  := -dynamiclib -Wl,-undefined,dynamic_lookup -Wl,-headerpad_max_install_names
  MATLAB_LIB_EXT  := dylib
  MATLAB_ARCH_DIR ?= maca64
  MEX_EXT ?= mex$(MATLAB_ARCH_DIR)
  MEX_CC ?= $(CXX)
  MEX_CXX ?= $(CXX)
  MATLAB_LIBOMP_REL := bin/$(MATLAB_ARCH_DIR)/libomp.$(MATLAB_LIB_EXT)
endif

OPENBLAS_SO := $(OPENBLAS_LIBDIR)/libopenblas.$(SHLIB_EXT)

PVFMM_A     := $(PVFMM_PREFIX)/lib/libpvfmm.a
FFTW_A_3    := $(FFTW_LIBDIR)/libfftw3.a
FFTW_A_3F   := $(FFTW_LIBDIR)/libfftw3f.a
OPENBLAS_A  := $(OPENBLAS_LIBDIR)/libopenblas.a
STATIC_LIBS := $(PVFMM_A) $(FFTW_A_3) $(FFTW_A_3F) $(OPENBLAS_A)

BACKEND := openblas

ifeq ($(BACKEND),accelerate)
  CPPFLAGS += -DUSE_ACCELERATE -DACCELERATE_NEW_LAPACK
  LDFLAGS  += -framework Accelerate
else ifeq ($(BACKEND),openblas)
#   # Prefer Homebrew's single-thread OpenBLAS if present; override with OPENBLAS_PREFIX if needed.
#   OPENBLAS_PREFIX ?= $(strip $(shell brew --prefix openblas-singlethread 2>/dev/null)) # /opt/homebrew/opt/openblas-singlethread
#   ifeq ($(strip $(OPENBLAS_PREFIX)),)
#     OPENBLAS_PREFIX := $(strip $(shell brew --prefix openblas 2>/dev/null)) # /opt/homebrew/opt/openblas
#   endif
#   OPENBLAS_PREFIX := $(strip $(OPENBLAS_PREFIX))
# #   OPENBLAS_PREFIX := /opt/homebrew/opt/openblas
#   # Always bring gfortran runtime when using OpenBLAS.
#   GFORTRAN_DIR ?= $(strip $(shell dirname $$($(GFORTRAN) -print-file-name=libgfortran.dylib 2>/dev/null))) # /opt/homebrew/opt/gcc/lib/gcc/current
#   GFORTRAN_DIR := $(strip $(GFORTRAN_DIR))
#   #   
#   CPPFLAGS += -I$(OPENBLAS_PREFIX)/include
#   LDFLAGS  += -L$(OPENBLAS_PREFIX)/lib -lopenblas \
#               -L$(GFORTRAN_DIR) -lgfortran -Wl,-rpath,$(GFORTRAN_DIR) \
#               -lquadmath
# # 	LDFLAGS  += $(OPENBLAS_PREFIX)/lib/libopenblas.a \
# #               -L$(GFORTRAN_DIR) -lgfortran -Wl,-rpath,$(GFORTRAN_DIR) \
# #               -lquadmath
  GFORTRAN_DIR := $(strip $(shell dirname $$($(GFORTRAN) -print-file-name=$(GFORTRAN_LIB) 2>/dev/null)))
  CPPFLAGS += -I$(OPENBLAS_PREFIX)/include
  CPPFLAGS += -I$(FFTW_PREFIX)/include
  LDFLAGS  := -L$(GFORTRAN_DIR) -lgfortran -lquadmath \
              -Wl,-rpath,$(GFORTRAN_DIR)
else
  $(error Unknown BACKEND '$(BACKEND)'; use accelerate or openblas)
endif

FFTW_STATIC_LDFLAGS := $(FFTW_A_3) $(FFTW_A_3F)
PVFMM_STATIC_LDFLAGS := $(PVFMM_A) -Wl,-rpath,$(PVFMM_PREFIX)/lib

.PHONY: help all lib mex pvfmm_octave_mex clean distclean

help:
	@echo "Available tasks:"
	@echo "  make                                      # runs lib, then mex, in order"
	@echo "  make lib                                  # build and install pvfmm into $(PVFMM_PREFIX)"
	@echo "  make mex [BACKEND=accelerate|openblas]     # build MATLAB MEX via mwrap/mex"
	@echo "  make pvfmm_octave_mex                     # build Octave MEX using g++"
	@echo "  make clean                                # remove generated binaries and gateways"
	@echo "  make distclean                            # clean + wipe pvfmm's build-omp/_install (forces a real rebuild)"
	@echo ""

all: lib mex

ifeq ($(UNAME_S),Linux)
  # Set MPI paths explicitly for Linux; adjust if mpi.h lives elsewhere.
  # Check with: `mpicc -show` or `rg --files /usr/include /usr/lib/x86_64-linux-gnu/openmpi/include | rg "mpi.h$"`.
  MPI_PREFIX := /usr/lib/x86_64-linux-gnu/openmpi
  MPI_LIBS   := -lmpi -lmpi_cxx
  OCTAVE_MPI_LIBS := -lmpi -lmpi_cxx
else
  MPI_PREFIX := /opt/homebrew/opt/open-mpi
  # MPI_LIBS   := -lmpi -lmpi_cxx
  MPI_LIBS   := -lmpi
  OCTAVE_MPI_LIBS := -lmpi
endif
MPI_INC := $(MPI_PREFIX)/include

PVFMM_CMAKE_ARGS = -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$(PVFMM_PREFIX) \
  -DCMAKE_C_COMPILER=$(CMAKE_C_COMPILER) \
  -DCMAKE_CXX_COMPILER=$(CMAKE_CXX_COMPILER) \
  -DCMAKE_C_FLAGS="-fPIC" \
  -DCMAKE_CXX_FLAGS="-fPIC" \
  -DCMAKE_PREFIX_PATH="$(OPENBLAS_PREFIX);$(FFTW_PREFIX)" \
  -DBLAS_LIBRARIES=$(OPENBLAS_SO) \
  -DLAPACK_LIBRARIES=$(OPENBLAS_SO) \
  -DMPI_HOME=$(MPI_PREFIX)

lib:
	$(CMAKE) -S $(PVFMM_SRC_DIR) -B $(PVFMM_BUILD_DIR) $(PVFMM_CMAKE_ARGS)
	$(CMAKE) --build $(PVFMM_BUILD_DIR) -j
	$(CMAKE) --build $(PVFMM_BUILD_DIR) --target install -j

# MWrap / MEX settings
MWRAP_MAC ?= /Users/hzhu/mwrap/mwrap
MWRAP_LINUX ?= ~/mwrap/mwrap
ifeq ($(UNAME_S),Linux)
  MWRAP ?= $(MWRAP_LINUX)
else
  MWRAP ?= $(MWRAP_MAC)
endif
MATLAB_ROOT ?= $(shell ls -d /Applications/MATLAB_R*.app /usr/local/MATLAB/R* /opt/MATLAB/R* 2>/dev/null | sort | tail -n1)
MEX ?= $(MATLAB_ROOT)/bin/mex
MEXFLAGS ?= -largeArrayDims
MEX_CXXSTD ?= -std=c++17
MWRAP_MEX_NAME := pvfmm
MWRAP_GATEWAY := $(MWRAP_MEX_NAME).c
MWRAP_MEX_OUTPUT := $(MWRAP_MEX_NAME).$(MEX_EXT)
MATLAB_INCLUDE := $(MATLAB_ROOT)/extern/include
SIMULINK_INCLUDE := $(MATLAB_ROOT)/simulink/include

MPI_CPPFLAGS := -I$(MPI_INC)
MPI_RPATH    := -Wl,$(comma)-rpath,$(comma)$(MPI_PREFIX)/lib
MPI_LDFLAGS  := -L$(MPI_PREFIX)/lib $(MPI_LIBS) $(MPI_RPATH)
OCTAVE_MPI_LDFLAGS := -L$(MPI_PREFIX)/lib $(OCTAVE_MPI_LIBS) $(MPI_RPATH)
PVFMM_CPPFLAGS := -I$(PVFMM_PREFIX)/include/pvfmm -I$(PVFMM_PREFIX)/examples/include

WRAP_SRC := pvfmm_wrap.cpp
IMPL_SRC :=
MWRAP_OBJ := $(MWRAP_MEX_NAME)_gw.o
PVFMM_WRAP_OBJ := pvfmm_wrap.o
PVFMM_C_OBJ := pvfmm_c.o

# Link the SAME OpenMP runtime MATLAB loads (install-name @rpath/libomp.dylib), NOT
# GNU libgomp -- otherwise two OpenMP runtimes in one process => at >1 thread the mex
# HANGS (the libgomp/libomp conflict; cf. make.inc.fmm3d). Thread count is then set by
# MATLAB's maxNumCompThreads(N), not OMP_NUM_THREADS (LLVM libomp ignores that env var).
MATLAB_LIBOMP := $(MATLAB_ROOT)/$(MATLAB_LIBOMP_REL)

.PHONY: mex
# Two-step build: COMPILE the gateway with -fopenmp so the #pragma omp regions in
# pvfmm_wrap.h emit GOMP_* calls; then LINK WITHOUT -fopenmp (drops the implicit
# -lgomp) and resolve those GOMP_* symbols from MATLAB's libomp instead -> exactly ONE
# OpenMP runtime. pvfmm_wrap.h is injected via -include (single TU, no separate .cpp).
mex: $(MWRAP_GATEWAY)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(MPI_CPPFLAGS) $(MEX_CXXSTD) -DNDEBUG \
	  $(PVFMM_CPPFLAGS) -I$(MATLAB_INCLUDE) -I$(SIMULINK_INCLUDE) \
	  -c $(MWRAP_GATEWAY) -o $(MWRAP_OBJ)
	$(CXX) $(filter-out -fopenmp,$(CXXFLAGS)) $(MEX_LINK_FLAGS) \
	  -o $(MWRAP_MEX_OUTPUT) \
	  $(MWRAP_OBJ) \
	  $(STATIC_LIBS) \
	  $(LDFLAGS) \
	  $(MPI_LDFLAGS) \
	  $(MATLAB_LIBOMP) \
	  -Wl,-rpath,$(MATLAB_ROOT)/bin/$(MATLAB_ARCH_DIR) \
	  -Wl,-rpath,$(MATLAB_ROOT)/extern/bin/$(MATLAB_ARCH_DIR) \
	  $(MATLAB_ROOT)/extern/bin/$(MATLAB_ARCH_DIR)/libMatlabEngine.$(MATLAB_LIB_EXT) \
	  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH_DIR)/libMatlabDataArray.$(MATLAB_LIB_EXT) \
	  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH_DIR)/libmex.$(MATLAB_LIB_EXT) \
	  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH_DIR)/libmx.$(MATLAB_LIB_EXT)

OCTAVE_CONFIG ?= octave-config
OCTAVE_INC ?= $(shell $(OCTAVE_CONFIG) -p OCTINCLUDEDIR 2>/dev/null)
OCTAVE_LIBDIR ?= $(shell $(OCTAVE_CONFIG) -p OCTLIBDIR 2>/dev/null)
OCTAVE_OCTFILEDIR ?= $(shell $(OCTAVE_CONFIG) -p OCTFILEDIR 2>/dev/null)
OCTAVE_LINK_OPTS ?= $(shell mkoctfile -p OCTAVE_LINK_OPTS 2>/dev/null)
OCTAVE_MEX_LIBS ?= $(shell mkoctfile -p LIBOCTMEX 2>/dev/null) $(shell mkoctfile -p OCTAVE_LIBS 2>/dev/null)
OCTAVE_MEX_EXT ?= mex
OCTAVE_CPPFLAGS := -I$(OCTAVE_INC) -DOCTAVE_MEX_FILE
OCTAVE_LDFLAGS := -L$(OCTAVE_LIBDIR) $(OCTAVE_MEX_LIBS) -Wl,-rpath,$(OCTAVE_LIBDIR)
ifneq ($(strip $(OCTAVE_OCTFILEDIR)),)
  OCTAVE_LDFLAGS += -Wl,-rpath,$(OCTAVE_OCTFILEDIR)
endif
ifeq ($(UNAME_S),Darwin)
  OCTAVE_LINK_OPTS := -bundle -undefined dynamic_lookup -Wl,-flat_namespace -Wl,-all_load
endif

.PHONY: pvfmm_octave_mex
pvfmm_octave_mex: $(MWRAP_GATEWAY) $(PVFMM_C_OBJ)
	@echo "Linking Octave MEX via mkoctfile with g++..."
	CC=$(CC) CXX=$(CXX) mkoctfile --mex -v \
	  -o $(MWRAP_MEX_NAME).mex \
	  $(PVFMM_C_OBJ) \
	  $(PVFMM_A) $(FFTW_A_3) $(FFTW_A_3F) $(OPENBLAS_A) \
	  $(OCTAVE_MPI_LDFLAGS) \
	  -lstdc++ \
	  -L$(GFORTRAN_DIR) -lgfortran -lquadmath -Wl,-rpath,$(GFORTRAN_DIR) \
	  -L$(MPI_PREFIX)/lib $(OCTAVE_MPI_LIBS) -Wl,-rpath,$(MPI_PREFIX)/lib

# pvfmm_mex: $(MWRAP_GATEWAY) $(WRAP_SRC)
# 	$(MEX) $(MEXFLAGS) \
# 		CC="$(MEX_CC)" \
# 		CXX="$(MEX_CXX)" \
# 		CFLAGS="$$CFLAGS $(CPPFLAGS)" \
# 		CXXFLAGS="$$CXXFLAGS $(MEX_CXXSTD) $(CPPFLAGS)" \
# 		LDFLAGS="$$LDFLAGS $(LDFLAGS) -Wl,-bundle" \
# 		-output $(MWRAP_MEX_NAME) \
# 		$(WRAP_SRC) $(IMPL_SRC)

$(MWRAP_GATEWAY): pvfmm.mw
	$(MWRAP) -mex $(MWRAP_MEX_NAME) -mb -list pvfmm.mw
	$(MWRAP) -mex $(MWRAP_MEX_NAME) -c $(MWRAP_GATEWAY) pvfmm.mw

# Compile the gateway as C++ (NOT gcc -std=c99): pvfmm_wrap.h is injected via
# -include and carries the C++ implementation, so the gateway TU must be C++.
$(PVFMM_C_OBJ): pvfmm.c
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(MPI_CPPFLAGS) $(MEX_CXXSTD) $(PVFMM_CPPFLAGS) $(OCTAVE_CPPFLAGS) -c pvfmm.c -o $@

clean:
	rm -f *.o $(MWRAP_MEX_NAME).$(OCTAVE_MEX_EXT) $(MWRAP_MEX_OUTPUT) $(MWRAP_GATEWAY) $(MWRAP_MEX_NAME)_gw.h

# Separate from `clean` on purpose: this wipes PVFMM's OWN build tree (build-omp,
# ~35 MB of CMakeFiles + .o) and its installed headers/lib (_install, ~8 MB) so the
# next `make lib` is a genuine from-scratch CMake configure + compile, not an
# incremental one that CMake may skip as up to date. Does NOT touch the pvfmm
# submodule's tracked source (git-safe; only removes untracked build output).
distclean: clean
	rm -rf $(PVFMM_BUILD_DIR) $(PVFMM_PREFIX)
