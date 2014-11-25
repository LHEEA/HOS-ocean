LOCDIR =.

SRCDIR =$(LOCDIR)/
SRCDIR2=$(LOCDIR)/../utilities/
OBJDIR =$(LOCDIR)/obj/
LIBDIR=/usr/local/lib/
BINDIR=$(LOCDIR)/bin/
## ifort compiler
#FC=ifort
##
#OPTFLAGS= -O3 # -fp-model precise this option to keep exactly same fp description
#DBFLAGS= -O0 -traceback -check all -warn all
#FLAGMOD1= -module $(OBJDIR) #Flag for writing modules in $(OBJ)
#FLAGMOD2= -module $(OBJDIR) #Flag for reading modules in $(OBJ)
#
# gfortran compiler (optimized for Mac Pro with intel corei7: avx not working on Mac)
FC=gfortran
#
OPTFLAGS= -O3 -march=corei7 -msse2 -funroll-loops -fno-protect-parens -ffast-math
DBFLAGS = -O0 -Wline-truncation -fpack-derived -fbacktrace -ffpe-summary=all -fimplicit-none -fcheck=all -Wall -Wtabs -Wextra -Wunderflow  -Wno-zerotrip -finit-integer=inf -finit-real=nan -ffpe-trap=zero,overflow
FLAGMOD1= -J $(OBJDIR) #Flag for writing modules in $(OBJ)
FLAGMOD2= -I $(OBJDIR) #Flag for reading modules in $(OBJ)

#define all modules
nomcode=Post_processing

module=\
	variables_post_process\
	input_post_process\
	read_files\
	analysis_wavefield\
	output_post_process\
	reconstruction

module2=\
	type\
	fourier_r2c_FFTW3

subroutine=\

# define the list of source files.
SRCS = $(addprefix $(SRCDIR), $(addsuffix .f90, $(nomcode))) $(addprefix $(SRCDIR), $(addsuffix .f90, $(module))) $(addprefix $(SRCDIR2), $(addsuffix .f90, $(module2)))

build:Release

# top-level rule, to compile everything.
all:clean depend Release

# Targets for compilation

Release: FFLAGS1 = $(OPTFLAGS) $(FLAGMOD1)
Release: FFLAGS2 = $(OPTFLAGS) $(FLAGMOD2)
Release: $(BINDIR)$(nomcode) $(SRCS)

Debug: FFLAGS1 = $(DBFLAGS) $(FLAGMOD1)
Debug: FFLAGS2 = $(DBFLAGS) $(FLAGMOD2)
Debug: $(BINDIR)$(nomcode) $(SRCS)

# now add a line to include the dependency list.
# after depend so that make depend is working without $(OBJDIR).depend
include $(OBJDIR).depend

Objlink = $(addprefix $(OBJDIR), $(addsuffix .o, $(module))) $(addprefix $(OBJDIR), $(addsuffix .o, $(module2))) $(addprefix $(OBJDIR), $(addsuffix .o, $(subroutine)))

# rule to link the program
$(BINDIR)$(nomcode): $(SRCDIR)$(nomcode:%=%.f90) $(Objlink) $(OBJDIR).depend
	$(FC) $(FFLAGS2) -o $(BINDIR)$(nomcode) $(SRCDIR)$(nomcode:%=%.f90) $(Objlink) $(LIBDIR)libfftw3.a $(LIBDIR)librefblas.a $(LIBDIR)liblapack.a

# now comes a meta-rule for compiling any "f90" source file.
$(OBJDIR)%.o: $(SRCDIR)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR2)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

clean:
	rm -f $(OBJDIR)*
	rm -f $(OBJDIR).depend
	rm -f $(BINDIR)$(nomcode)
	@echo "all cleaned up!"

# rule for building dependency lists, and writing them to a file
# named ".depend".
depend $(OBJDIR).depend:
	@if test -d bin; then echo "bin/ exists"; else mkdir bin; fi
	@if test -d obj; then echo "obj/ exists"; else mkdir obj; fi
	rm -f $(OBJDIR).depend
	makedepf90 -b $(OBJDIR) $(SRCS) > $(OBJDIR).depend

.PHONY: clean depend

