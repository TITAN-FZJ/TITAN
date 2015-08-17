SHELL = /bin/sh
####################################################
####################################################
#             SHE makefile
####################################################
#
# to compile use the following options:
#
# make          		       -  default compilation using mpif90
#
# make PLATFORM=(location) -  uses specific compiler options
# Configured locations:
# uff, iff, osx, juropa, juropatest, jureca, juqueen
#
# make DEBUG=debug         -  uses debug flags
#
# make PARALLEL=omp        -  uses openmp
#
# make FILE=(filename.exe) -  creates executable (filename.exe)
#
# make clean    			     -  removes *.o *.mod *__genmod* *.exe files
#
# make cleanall  			     -  also remove dependency files .dep
#
####################################################
# Filename and folders configuration               #
####################################################
SRCDIR = ./source
BINDIR = ./$(LATTICE)/$(SYSTEM)
OBJDIR = ./$(LATTICE)/$(SYSTEM)/build
ifdef FILE
	FILENAME = $(BINDIR)/$(FILE)
else
	FILENAME = $(BINDIR)/main.exe
endif

####################################################
# Suffixes used in this makefile                   #
####################################################
.SUFFIXES:
.SUFFIXES: .F90 .o .mod .dep

####################################################
# Source files                                     #
####################################################
SRC =$(wildcard $(SRCDIR)/*.F90)

####################################################
# Objects to compile                               #
####################################################
OBJ = $(addprefix $(OBJDIR)/,$(notdir $(SRC:.F90=.o)))

####################################################
# Dependency files                                #
####################################################
DEP = $(OBJ:.o=.dep)

#=======================================================================
#============================ DEFAULT VALUES ===========================
#=======================================================================

####################################################
#  Compiler                                        #
####################################################
FC = mpif90

####################################################
#  Libraries                                       #
####################################################
LLIBS =-mkl -shared-intel -lnag

####################################################
#  Flags                                           #
####################################################
FFLAGS =-O3 -xSSE4.2

####################################################
#  Preprocessor                                    #
####################################################
CPP = -fpp

####################################################
#  Debugger Flags                                  #
####################################################
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif

####################################################
#  Module and Include folder                       #
####################################################
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/

####################################################
#  Parallelization                                 #
####################################################
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif

#=======================================================================
#========================== LOCATION SPECIFIC ==========================
#=======================================================================

####################################################
#                       UFF                        #
####################################################
ifeq ($(PLATFORM),uff)
# Compiler
FC = mpif90
# Preprocessor
CPP = -fpp -D _UFF
# Libraries
LLIBS =-mkl -static-intel -L$(HOME)/lib -lnag
#Flags
FFLAGS =-O3 -xSSE4.2
#Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif
####################################################
#                       OSX                        #
####################################################
ifeq ($(PLATFORM),osx)
# Compiler
FC = /usr/local/impi/bin/mpif90
# Preprocessor
CPP = -fpp
# Libraries
LLIBS =-mkl -shared-intel -L$(HOME)/lib -lkibe
#Flags
FFLAGS =-O3 -xSSE4.2
#Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif
####################################################
#                       IFF                        #
####################################################
ifeq ($(PLATFORM),iff)
# Compiler
FC = mpiifort
# Preprocessor
CPP = -fpp
# Libraries
LLIBS =-mkl -shared-intel -L$(HOME)/lib -lnag
# Flags
FFLAGS =-O3 -xSSE4.2
# Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif
####################################################
#                     JUROPA                       #
####################################################
ifeq ($(PLATFORM),juropa)
# Compiler
FC = mpif90
# Preprocessor
CPP = -fpp
# Libraries
LLIBS =-mkl -shared-intel -L$(HOME)/lib -lkibe
#Flags
FFLAGS =-O3 -xSSE4.2
# Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif
####################################################
#                   JUROPATEST                     #
####################################################
ifeq ($(PLATFORM),juropatest)
# Compiler
FC = mpif90
# Preprocessor
CPP = -fpp
# Libraries
LLIBS =-mkl -shared-intel -L$(HOME)/lib -lkibe
#Flags
FFLAGS =-O3 -xSSE4.2
#Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif
####################################################
#                     JURECA                       #
####################################################
ifeq ($(PLATFORM),jureca)
# Compiler
FC = mpiifort
# Preprocessor
CPP = -fpp
# Libraries
LLIBS =-mkl -shared-intel -L$(HOME)/lib -lkibe
#Flags
FFLAGS =-O3 -xSSE4.2
#Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-CB -check all -check uninit -ftrapuv -debug all -traceback -g -warn all -O0
endif
FFLAGS += -module $(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-openmp
FFLAGS +=-openmp
endif
endif
####################################################
#                    JUQUEEN                       #
####################################################
ifeq ($(PLATFORM),juqueen)
# Compiler
FC = mpif90
# Preprocessor
CPP = -WF,-D_JUQUEEN
# Libraries
LLIBS =$(ND)/lib/libnag_essl.a -L/bgsys/local/lib -lesslbg
#LLIBS =-L/opt/ibmmath/essl/5.1/lib64 -lesslbg -L$(HOME)/lib -lkibe
#Flags
FFLAGS =-O3 -qstrict -qnoescape -qnosave
#Debugger
ifeq ($(DEBUG),debug)
FFLAGS =-C -g -O0 -qflttrap -qinitauto=7FF7FFFF -qnoescape
endif
FFLAGS += -qmoddir=$(OBJDIR)/ -I$(OBJDIR)/
# Parallelization
ifeq ($(PARALLEL),omp)
LLIBS +=-qsmp=omp -qthreaded
FFLAGS +=-qsmp=omp -qthreaded
endif
endif
####################################################
# Linking                                          #
####################################################
all: $(FILENAME)

$(FILENAME): $(OBJ)
	@echo Creating executable $@ $(and $(strip $(DEBUG) $(PARALLEL)), with $(strip $(DEBUG) $(PARALLEL)))
	@$(FC) $^ -o $@ $(LLIBS)

$(OBJDIR)/%.o $(OBJDIR)/%.mod: $(SRCDIR)/%.F90
	@echo Compiling file $< to $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.F90,%.o,$<))) $(and $(strip $(DEBUG) $(PARALLEL)), with $(strip $(DEBUG) $(PARALLEL)))
	@$(FC) $(FFLAGS) $(CPP) -c $< -o $(addprefix $(OBJDIR)/,$(notdir $(patsubst %.F90,%.o,$<)))

####################################################
# Dependencies                                     #
####################################################
$(OBJDIR)/%.dep: $(SRCDIR)/%.F90 f90_mod_deps.py
	@echo Building dependency of file $<
	@./f90_mod_deps.py -o $@ -d "(mod_.*)" -D "$(OBJDIR)/\1.mod" -m "(.*)" -M "$(OBJDIR)/\1.mod" -O "$(OBJDIR)/\1.o" $<

ifeq ($(filter $(MAKECMDGOALS),clean cleanall),)
-include $(DEP)
endif

####################################################
# Clean up                                         #
####################################################
clean:
	@echo Deleting executable $(FILENAME) and cleaning folder $(OBJDIR)/
	@rm -f $(OBJDIR)/*
	@rm -f $(FILENAME).ld* $(FILENAME)

# cleanall:
# 	rm -f $(SRCDIR)/*.dep
# 	rm -f $(LATDIR)/*.dep
# 	rm -f $(SYSDIR)/*.dep
# 	rm -f $(OBJDIR)/*
# 	rm -f $(FILENAME).ld* $(FILENAME)