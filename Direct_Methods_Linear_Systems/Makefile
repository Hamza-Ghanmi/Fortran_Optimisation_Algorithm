#
# This is the Makefile for the model
#

# Program Name
PROG = main.out

# Source Folder name
VPATH=src
# Object Folder
OBJDIR=objs
# Module Folder
MODDIR=mods
# Executable Folder
EXECDIR=exec

# Compiler and Flags
FC = /usr/bin/gfortran
FFLAGS = -c -O3 -mcmodel=medium -g -I/usr/include
FLINK = -O3 -mcmodel=medium -g -L/usr/lib -lnetcdff
LINKER = $(FC) -o

# Object files
OBJS = main.o determinant.o cramer.o substitution.o Gaussian_elimination.o LU.o Cholesky.o

model: $(PROG)

# Creates the model
$(PROG): $(OBJS)
	@echo "--------------------------------------"
	@echo "Creating the executable for the model"
	@echo "--------------------------------------"
	$(LINKER) $(PROG) $(OBJS) $(FLINK)
	mv *.o $(OBJDIR)
	mv *.mod $(MODDIR)
	mv *.out $(EXECDIR)

%.o: %.f95
	@echo "--------------------------------------"
	@echo "Compiling the file $<"
	@echo "--------------------------------------"
	$(FC) $(FFLAGS) $<

# Cleans up everything
clean:
	@echo "--------------------------------------"
	@echo "Cleaning everything up in model"
	@echo "--------------------------------------"
	rm -f *~ *.nc plot*.png *.out *.o *.mod
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*~
	rm -f $(MODDIR)/*.mod $(MODDIR)/*~
	rm -f $(EXECDIR)/*~ $(EXECDIR)/*.out $(EXECDIR)/*.nc
	rm -f $(VPATH)/*~

main.o : main.f95 determinant.o cramer.o substitution.o
determinant.o : determinant.f95
cramer.o : cramer.f95
substitution.o : substitution.f95
aussian_elimination.o : aussian_elimination.f95
LU.o : LU.f95
Colesky.o : Cholesky.f95
