## subroutine and function source files
SRC_FILES = 
## module files
MDL_FILES = mkl_spblas.f90 globaldata.f90 tools.f90  io.f90 macroscopic.f90 microscopic_flux.f90 microscopicSn.f90
#MDL_FILES = mkl_spblas.f90 globaldata.f90 tools.f90  io.f90 macroscopic.f90 microscopic_part.f90 microscopicSn.f90
OBJS = ${SRC_FILES:.f90=.o}
MOD = ${MDL_FILES:.f90=.o}
## main file
BIN = main

# compiler
F90 = mpiifort
FFLAGS = -qmkl -qopenmp

##顺序非常重要，决定了依赖关系
OBJECTS =$(MOD) $(OBJS) $(BIN).o
SOURCES = $(SRC_FILES) $(MDL_FILES) $(BIN).f90

## final target: note the "tab" in the command line
$(BIN): $(OBJECTS)
	$(F90) $(FFLAGS) -o $(BIN) $(OBJECTS)
$(OBJS): $(SRC_FILES)
	$(F90) $(FFLAGS) -c $(SRC_FILES)
$(MOD): $(MDL_FILES)
	$(F90) $(FFLAGS) -c $(MDL_FILES)
                                                                                    
$(BIN).o: $(SOURCES)
	$(F90) $(FFLAGS) -c $(BIN).f90
clean:
	rm -f *.plt *.o *.mod *.out $(BIN)
