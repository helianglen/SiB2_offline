#
# Makefile for SIB2 model pre-process
#
# Last revision: 20151030 by nelsonvn
#
# Usage:    make COMPILER=opt_comp CPU=opt_cpu BUILD=opt_build
#
# where 
# opt_comp  = [gnu|intel|amd|pgi]
# opt_cpu   = [amd|intel]
# opt_build = [debug|fast]
#

# Default options
COMPILER = gnu
CPU = intel
BUILD = debug

include paths.mk
include include.mk.$(COMPILER)

# This line replaces the OBJDIR and MODDIR default
OBJDIR = ./obj
MODDIR = ./obj

MODEL_SRC = \
	constants.f90 \
	basic_vars.f90 \
	cal2jul.f90 \
	common_inc.f90 \
	comsibc_h.f90 \
	pardif_h.f90 \
	sib2par_inc.f90 \
	sib2data.f90 \
	SIB2sub.f90 \
	derive_trans.f \
	sib2_offline.f90
#	sib2_sub.f90 

MODEL_OBJ = $(patsubst %.f90,%.o, $(patsubst %.f,%.o,$(MODEL_SRC)))

FC = $(F77)

%.o : %.f90
	$(FC) $(FFLAGS) -c $< -o $(OBJDIR)/$(@F)

%.o : %.f
	$(FC) $(FFLAGS) -c $< -o $(OBJDIR)/$(@F)

all:	$(MODEL_OBJ)
	$(FC) $(FFLAGS) -o sib2_offline.out $(OBJDIR)/*.o

clean:
	rm -f $(MODDIR)/* $(OBJDIR)/*

show_mod:
#	The following command show the modules source list
	@echo $(MODULES_SRC)

show_src:
#	The following command show the model source code list
	@echo $(MODEL_SRC)

tar:
	rm -f sib2_offline.tar.bz2
	tar -jcf sib2_offline.tar.bz2 $(MODEL_SRC) data2 *.par *.sh \
	Makefile paths.mk include.*
