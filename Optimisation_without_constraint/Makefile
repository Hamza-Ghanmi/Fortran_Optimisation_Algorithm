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
OBJS = main.o gdfs.o gdos.o gdcs.o

model: $(PROG)

# Creates the model
$(PROG): $(OBJS)
	@echo "--------------------------------------"
	@echo "Creating the executable for the model"
	@echo "--------------------------------------"
	$(LINKER) $(PROG) $(OBJS) $(FLINK)
	mv *.o $(OBJDIR)
	#mv *.mod $(MODDIR)
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

main.o : main.f95 gdfs.o
gdfs.o : gdfs.f95
gdos.o : gdos.f95
gdcs.o : gdcs.f95
