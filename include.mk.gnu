#
# include.mk.gnu
#

# Fortran compiler
F77 = gfortran

# Flags for AMD processors
ifeq ($(CPU),amd)
   ifeq ($(BUILD), debug)
#      FFLAGS = -ffixed-line-length-none -J$(MODDIR) -march=bdver1 -O2 -Warray-bounds -fbacktrace -std=f95
      FFLAGS = -ffixed-line-length-none -J$(MODDIR) -march=bdver1 -O2 \
	       -fbacktrace \
      	       -Warray-bounds \
	       -Wunused \
      	       -Wuninitialized \
	       -Wconversion \
	       -finit-local-zero \
	       -fdefault-real-8
   endif

   ifeq ($(BUILD),fast)
      FFLAGS = -ffixed-line-length-none -J$(MODDIR) -march=bdver1 -O3 \
	       -mprefer-avx128 \
	       -funroll-loops \
               -minline-all-stringops \
	       -mieee-fp \
	       -ftree-vectorize \
	       -fdefault-real-8
   endif
endif

# Flags for INTEL processors
ifeq ($(CPU),intel)
   ifeq ($(BUILD),debug)
#      FFLAGS = -ffixed-line-length-none -J$(MODDIR) -march=native -O2 -Warray-bounds -frange-check -fbacktrace
#      FFLAGS = -ffixed-line-length-none -J$(MODDIR) -g -O0 -finit-integer=-99999 -finit-real=nan -fdefault-real-8
      FFLAGS = -ffixed-line-length-none -J$(MODDIR) -g -O0 -fdefault-real-8 \
   	       -fbacktrace \
	       -Wunused \
      	       -Wuninitialized \
	       -Wunused-dummy-argument \
      	       -Warray-bounds \
	       -Wconversion \
	       -fcheck=all
#	       -fdefault-real-8 \
#	       -finit-local-zero
   endif

   ifeq ($(BUILD),fast)
      FFLAGS = -ffixed-line-length-none -J$(MODDIR) -march=native -O3 \
      	       -ffast-math \
	       -funroll-loops \
               -fipa-cp \
	       -fdefault-real-8
   endif
endif
