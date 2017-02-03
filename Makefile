
# ******************** JK6 MAKEFILE **************** 

#----------------------------------------------------------------------------
#                         USER DEFINABLE OPTIONS
#----------------------------------------------------------------------------

# Executable name
EXENAME = leporello

# Print relevant output to log file
LOGFILE = yes
LOGNAME = leporello.log

# Compiler ( gfortran, ifort )
FC = ifort

# Debugging options ( yes or no )
DEBUG =  yes 

# Optimization level
OPTLEVEL = 3

# linking FFTW 3.3
FFTW3 = no 

# OpenMP libraries
OPENMP = no 

# linking LAPACK and BLAS 
LAPACK = yes

# Compile with standard real 8 (see details about the flags for each compiler...)
REAL8 = yes

#----------------------------------------------------------------------------
#                      LAPACK AND BLAS FINE DETAILS
#----------------------------------------------------------------------------

# Intel compiler version ( used only if FC=ifort and LAPACK=yes )
# 2016-SEQ, 2016-MULTI      -   2016 version, sequential / multithreaded 
# 2013-SEQ, 2013-MULTI      -   2013 version, sequential / multithreaded 
# 11-SEQ,   11-MULTI        -   11.x version, sequential / multithreaded 
# 11-IA32                   -   11.x ia32 arch version, sequential  
# 10-SEQ,   10-MULTI        -   10.x version, sequential / multithreaded 
INTELVERS = 2016-SEQ

# gfortran lapack libraries
# GNU      - system default libraries 
# ATLAS    - atlas libraries
GLAPACK = GNU

#----------------------------------------------------------------------------
#                             STRIP ALL SPACES
#----------------------------------------------------------------------------

# Strip leading and trailing spaces from all variables.
FC := $(strip ${FC})
DEBUG := $(strip ${DEBUG})
OPTLEVEL := $(strip ${OPTLEVEL})
LAPACK := $(strip ${LAPACK})
INTELVERS := $(strip ${INTELVERS})
GLAPACK := $(strip ${GLAPACK})
LOGFILE := $(strip ${LOGFILE})
LOGNAME := $(strip ${LOGNAME})
FFTW3 := $(strip ${FFTW3})
OPENMP := $(strip ${OPENMP})


#----------------------------------------------------------------------------
#                       Compiler specific statements.
#----------------------------------------------------------------------------

ifeq (${FC},gfortran)

   # Optimization flag
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flag
   DEBUGFLG = -g -fbounds-check

   # LAPACK AND BLAS flags
   ifeq (${GLAPACK},GNU)
      # GNU Lapack and Blas flags 
      LAPACKFLG = -llapack -lblas
   endif
   ifeq (${GLAPACK},ATLAS)
      # ATLAS Lapack and Blas flags
      LAPACKFLG =  -L/usr/lib64/atlas/ -llapack -lf77blas -lcblas -latlas
   endif

   # FFTW3 flags
   FFTW3FLG = -L/usr/local/lib/ -lfftw3
   FFTW3COMPILE = -I/usr/local/include/ 

   # OPENMP flags
   OPENMPFLG = -fopenmp

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
	DATAFLG = -fdefault-real-8
   endif

   # Flag to specify the position of mod files
   MODULEFLG = -fintrinsic-modules-path

endif

ifeq (${FC},ifort)

   # Optimization flags
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flags
   DEBUGFLG  =  -g -traceback -fpe-all=0 -debug all -check all 

   # MKL flags
   ifeq (${INTELVERS},2016-SEQ)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -mkl=sequential
   endif
   ifeq (${INTELVERS},2016-MULTI)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -qopenmp -mkl=parallel
   endif
   ifeq (${INTELVERS},2013-SEQ)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -mkl=sequential
   endif
   ifeq (${INTELVERS},2013-MULTI)
      LAPACKFLG = -lpthread -lm
      LAPACKCOMPILE = -openmp -mkl=parallel 
   endif
   ifeq (${INTELVERS},11-SEQ)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
      LAPACKCOMPILE = -I$(MKLROOT)/include
   endif
   ifeq (${INTELVERS},11-IA32)
      LAPACKFLG = -L${MKLROOT}/lib/ia32 -lmkl_intel -lmkl_core -lmkl_sequential -lpthread -lm 
      LAPACKCOMPILE = -I${MKLROOT}/include
   endif
   ifeq (${INTELVERS},11-MULTI)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
      LAPACKCOMPILE = -openmp -I$(MKLROOT)/include
   endif
   ifeq (${INTELVERS},10-SEQ)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
      LAPACKCOMPILE = -I$(MKLROOT)/include
   endif
   ifeq (${INTELVERS},10-MULTI)
      LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
      LAPACKCOMPILE = -openmp -I$(MKLROOT)/include
   endif

   # FFTW3 flags
   FFTW3FLG = -L/data/local/fftw3/lib/ -lfftw3
   FFTW3COMPILE = -I/data/local/fftw3/include/ 

   # OPENMP flags
   OPENMPFLG = -qopenmp

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
	DATAFLG = -r8 -i4
   endif
 
   # Flag to specify the position of mod files
   MODULEFLG = -I

endif

ifeq (${FC},fermi)

   # Set name of the cross-compiler wrapper
   FC = bgxlf_r -qarch=qp -qtune=qp

   # Optimization flags
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flags
   DEBUGFLG  = -g -qfullpath -qcheck -qflttrap -qkeepparm

   # MKL flags
   LAPACKFLG = -L\${LAPACK_LIB} -llapack -L\${ESSL_LIB} -lesslbg -L\${BLAS_LIB} -lblas
   LAPACKCOMPILE =

   # FFTW3 flags
   ifeq (${REAL8},yes)
      FFTW3FLG = -L\${FFTW3_LIB} -lfftw3
      FFTW3COMPILE = -I${FFTW3_INC}
   else
      FFTW3FLG = -L\${FFTW3_LIB} -lfftw3f
      FFTW3COMPILE = -I\${FFTW3_INC}
   endif

   # OPENMP flags
   OPENMPFLG = -openmp

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
      DATAFLG = -r8 -i4
   endif

   # Flag to specify the position of mod files
   MODULEFLG = -I

endif

ifeq (${FC},hilbert)

   # Set name of the intel compiler
   FC = ifort

   # Optimization flags
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flags
   DEBUGFLG  =  -g -traceback -fpe-all=0 -debug all -check all

   # Intel 11 with threading
   LAPACKFLG = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm
   LAPACKCOMPILE = -openmp -I$(MKLROOT)/include

   # FFTW3 flags
   FFTW3FLG = -lfftw3

   # OPENMP flags
   OPENMPFLG = -openmp

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
        DATAFLG = -r8 -i4
   endif

   # Flag to specify the position of mod files
   MODULEFLG = -I

endif

ifeq (${FC},galileo)

   # Set name of the intel compiler
   FC = ifort
   
   # Optimization flags
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flags
   DEBUGFLG  =  -g -traceback -fpe-all=0 -debug all -check all
   
   # Intel mkl with threading
   LAPACKFLG = -lpthread -lm
   LAPACKCOMPILE = -mkl=sequential

   # FFTW3 flags
   FFTW3FLG = -L${FFTW_LIB} -lfftw3 
   FFTW3COMPILE = -I${FFTW_INC} 
  
   # OPENMP flags
   OPENMPFLG = -openmp

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
      DATAFLG = -r8 -i4
   endif
  
   # Flag to specify the position of mod files
   MODULEFLG = -I

endif



#----------------------------------------------------------------------------
#              Setup preprocessing options 
#----------------------------------------------------------------------------

PPDEFINE = 

# Define macro for writing current version (git tag) in the code
PPDEFINE += -DVERSIONTAG=\"$(shell git describe --tags)\"

# Setup preprocessing debug options
ifeq (${DEBUG}, yes)
   PPDEFINE += -DVERBOSE_OUTPUT
endif

# Write output to logfile
ifeq (${LOGFILE}, yes)
   PPDEFINE += -DLOG_FILE=\"${LOGNAME}\"
endif

# Preprocess with lapack calls
ifeq (${LAPACK}, yes)
   PPDEFINE += -DWITH_LAPACK
endif

# Preprocess using FFTW3 calls
ifeq (${FFTW3}, yes)
   PPDEFINE += -DWITH_FFTW3
endif

# Preprocess with OPENMP
ifeq (${OPENMP}, yes)
   PPDEFINE += -DWITH_OPENMP
endif


#----------------------------------------------------------------------------
#              Setup linking and compilation flags
#----------------------------------------------------------------------------

# initialize flags
COMPILEFLG = 
LINKFLG    = 
LIBFLG     = 
INCLUDEFLG =

# if debugging set the appropriate flags
ifeq (${DEBUG}, yes)
   COMPILEFLG += ${DEBUGFLG}
endif

# Set flags for defining standard variable kinds
COMPILEFLG += ${DATAFLG} 

# If lapack, add the linking options
ifeq (${LAPACK}, yes)
   LIBFLG += ${LAPACKFLG}
   LINKFLG += ${LAPACKCOMPILE}
endif

# If FFTW3, add the linking options
ifeq (${FFTW3}, yes)
   LIBFLG += ${FFTW3FLG}
   COMPILEFLG += ${FFTW3COMPILE}
endif

# If OPENMP, add the linking options
ifeq (${OPENMP}, yes)
   LIBFLG += ${OPENMPFLG}
   COMPILEFLG += ${OPENMPFLG}
endif


#----------------------------------------------------------------------------
#              Determine the optimization level to be used.
#----------------------------------------------------------------------------

# if debugging override input optimization level
ifeq (${DEBUG}, yes)
   OPTLEVEL = 0
endif

# Define optimization level
OPTFLAGS       = ${O0FLAGS}
ifeq (${OPTLEVEL},1)
  OPTFLAGS       = ${O1FLAGS}
endif
ifeq (${OPTLEVEL},2)
  OPTFLAGS       = ${O2FLAGS}
endif
ifeq (${OPTLEVEL},3)
  OPTFLAGS       = ${O3FLAGS}
endif

COMPILEFLG += ${OPTFLAGS}


#----------------------------------------------------------------------------
#                      List of directories
#----------------------------------------------------------------------------

SRCDIR  = Source
PPDIR   = PreProcessed
OBJDIR  = Objects
TESTDIR = Tests/Source
EXEDIR  = Executables

#----------------------------------------------------------------------------
#                      List of object files
#----------------------------------------------------------------------------

# Define list of object from the list of all f90 files in the directory
OBJSWITHMAIN =$(patsubst Source/%,Objects/%,$(patsubst %.f90,%.o,$(wildcard ${SRCDIR}/*.f90)))
OBJS =$(patsubst %/Main.o,,${OBJSWITHMAIN})


#----------------------------------------------------------------------------
#       Construct the compile, link, preprocess variables.
#----------------------------------------------------------------------------

# Compile command: ${COMPILE} <source>
COMPILE                 = ${FC} ${COMPILEFLG} ${MODULEFLG} ${OBJDIR} -c 

# Link command: ${LINK} <exe name> <objects> <libflags>
LINK                    = ${FC} ${LINKFLG} ${MODULEFLG} ${OBJDIR} -o

# Preprocess: ${PREPROCESS} <to preprocess> <preprocessed>
PREPROCESS              = cpp -P -traditional ${PPDEFINE}

# Build static library : ${AR} <libraryname>.a <objects>
AR 			= ar cr

#----------------------------------------------------------------------------
#                       START OF MAKE RULES
#----------------------------------------------------------------------------

# Link objects to the produce the executable file
${EXENAME} : ${SRCDIR}/Main.f90 ${OBJS} 
	${PREPROCESS} ${SRCDIR}/Main.f90 ${PPDIR}/Main.f90
	${COMPILE} ${PPDIR}/Main.f90 
	${LINK} ${EXEDIR}/$@ Main.o $(OBJS) ${LIBFLG}
	rm Main.o

# Make target to build all the object files and assemble them
all : ${OBJS}

# Make a target object file by preprocessing and compiling the fortran code
${OBJDIR}/%.o : ${SRCDIR}/%.f90
	${PREPROCESS} ${SRCDIR}/$*.f90 ${PPDIR}/$*.f90
	${COMPILE} ${PPDIR}/$*.f90 
	cp -p $*.o $(shell echo $* | tr A-Z a-z).mod ${OBJDIR}
	rm $*.o $(shell echo $* | tr A-Z a-z).mod

# Make target to build required directories
directories : 
	mkdir -p ${PPDIR} ${OBJDIR} ${EXEDIR}

# Make documentation with doxygen
doc :
	doxygen Documentation/Doxyfile

# Remove compiled objects and related stuff
clean :
	rm -fr ${OBJDIR}/* ${PPDIR}/* 

# Clean documentation
clean-doc :
	rm -fr Documentation/html
	rm -fr Documentation/latex

# objects with different rules
${OBJDIR}/vEleyRideal_3D.o : ${SRCDIR}/vEleyRideal_3D.f90
	${COMPILE} ${SRCDIR}/vEleyRideal_3D.f90 
	cp -p vEleyRideal_3D.o ${OBJDIR}
	rm vEleyRideal_3D.o

# --------------------------------------------------------------------------------------------
# ---------------------------------- START WITH DEPENDENCIES NOW -----------------------------
# --------------------------------------------------------------------------------------------

# Very basic files, which everything depends on
COMMONDEP = Makefile ${OBJDIR}/ErrorTrap.o  ${OBJDIR}/MyConsts.o ${OBJDIR}/MyLinearAlgebra.o  ${SRCDIR}/preprocessoptions.cpp

# Program to simulate the vibrational relaxation
#${OBJDIR}/VibrationalRelax.o : ${SRCDIR}/VibrationalRelax.f90 ${OBJDIR}/SharedData.o ${OBJDIR}/InputField.o                       \
#                                ${OBJDIR}/UnitConversion.o ${OBJDIR}/ClassicalEqMotion.o ${OBJDIR}/PotentialModule.o              \
#                                ${OBJDIR}/IndependentOscillatorsModel.o ${OBJDIR}/RandomNumberGenerator.o ${COMMONDEP}

# Program to simulate the thermal equilibrium dynamics
#${OBJDIR}/ThermalEquilibrium.o : ${SRCDIR}/ThermalEquilibrium.f90 ${OBJDIR}/SharedData.o ${OBJDIR}/InputField.o                    \
#                                  ${OBJDIR}/ClassicalEqMotion.o ${OBJDIR}/PotentialModule.o ${OBJDIR}/IndependentOscillatorsModel.o \
#                                  ${OBJDIR}/RandomNumberGenerator.o ${COMMONDEP}

# Program to simulate the scattering process
${OBJDIR}/ScatteringSimulation.o : ${SRCDIR}/ScatteringSimulation.f90 ${OBJDIR}/SharedData.o ${OBJDIR}/InputField.o ${OBJDIR}/UnitConversion.o  \
				${OBJDIR}/ClassicalEqMotion.o ${OBJDIR}/PotentialModule.o ${OBJDIR}/IndependentOscillatorsModel.o               \
				${OBJDIR}/RandomNumberGenerator.o ${OBJDIR}/PrintTools.o ${COMMONDEP}                                           \

# Program to simulate the vibrational relaxation
#${OBJDIR}/Harmonic1DModel.o : ${SRCDIR}/Harmonic1DModel.f90 ${OBJDIR}/SharedData.o ${OBJDIR}/InputField.o                         \
#                               ${OBJDIR}/ClassicalEqMotion.o ${OBJDIR}/IndependentOscillatorsModel.o                              \
#                               ${OBJDIR}/RandomNumberGenerator.o ${OBJDIR}/UnitConversion.o ${OBJDIR}/FFTWrapper.o ${COMMONDEP}

# Program to compute the properties of the PES
${OBJDIR}/PotentialAnalysis.o : ${SRCDIR}/PotentialAnalysis.f90  ${OBJDIR}/SharedData.o ${OBJDIR}/InputField.o                     \
				${OBJDIR}/UnitConversion.o ${OBJDIR}/PotentialModule.o  ${OBJDIR}/FiniteDifference.o               \
				${OBJDIR}/PrintTools.o  ${COMMONDEP}

# Program to compute the minimum energy path of the PES
#${OBJDIR}/MinimumEnergyPath.o : ${SRCDIR}/MinimumEnergyPath.f90 ${OBJDIR}/SharedData.o ${OBJDIR}/InputField.o                     \
#                                 ${OBJDIR}/UnitConversion.o ${OBJDIR}/PotentialModule.o ${OBJDIR}/IndependentOscillatorsModel.o   \
# 				${OBJDIR}/PrintTools.o ${COMMONDEP}

# Module containing the definitions of the independent oscillator model
${OBJDIR}/IndependentOscillatorsModel.o  : ${SRCDIR}/IndependentOscillatorsModel.f90 ${OBJDIR}/SplineInterpolator.o                \
				   ${OBJDIR}/RandomNumberGenerator.o ${OBJDIR}/UnitConversion.o ${OBJDIR}/FFTWrapper.o ${COMMONDEP}

# Module containing the integrator for the classical eq of motion
${OBJDIR}/ClassicalEqMotion.o  : ${SRCDIR}/ClassicalEqMotion.f90 ${OBJDIR}/RandomNumberGenerator.o ${OBJDIR}/FFTWrapper.o         \
 					${COMMONDEP}

# Module containing the potential energy surface
${OBJDIR}/PotentialModule.o  : ${SRCDIR}/PotentialModule.f90 ${OBJDIR}/vEleyRideal_3D.o ${OBJDIR}/RandomNumberGenerator.o          \
				${OBJDIR}/FiniteDifference.o ${COMMONDEP}

# Module containing the common data (v3)
${OBJDIR}/SharedData.o : ${SRCDIR}/SharedData.f90 ${OBJDIR}/IndependentOscillatorsModel.o ${COMMONDEP}

# Module containing the subroutine to write data in vtk format
#${OBJDIR}/PrintTools.o : ${SRCDIR}/PrintTools.f90 ${COMMONDEP}

# Module containing the random number generator
${OBJDIR}/RandomNumberGenerator.o : ${SRCDIR}/RandomNumberGenerator.f90 ${COMMONDEP}

# Input and output unit conversion
${OBJDIR}/UnitConversion.o : ${SRCDIR}/UnitConversion.f90 ${COMMONDEP}

# Set procedure for reading quasi-free format input file
${OBJDIR}/InputField.o : ${SRCDIR}/InputField.f90 ${COMMONDEP}

# Wrapper for FFTW 3.3
${OBJDIR}/FFTWrapper.o : ${SRCDIR}/FFTWrapper.f90 ${COMMONDEP}

# Compute derivatives with finite difference
${OBJDIR}/FiniteDifference.o : ${SRCDIR}/FiniteDifference.f90 ${OBJDIR}/RandomNumberGenerator.o ${COMMONDEP}

# Module containing the spline interpolation subroutines
${OBJDIR}/SplineInterpolator.o : ${SRCDIR}/SplineInterpolator.f90 ${OBJDIR}/NRUtility.o Makefile

# Linear algebra module
${OBJDIR}/MyLinearAlgebra.o : ${SRCDIR}/MyLinearAlgebra.f90 ${OBJDIR}/ErrorTrap.o ${OBJDIR}/MyConsts.o              \
				${OBJDIR}/NRUtility.o Makefile

# Utility subroutines and functions of NR
${OBJDIR}/NRUtility.o : ${SRCDIR}/NRUtility.f90 Makefile

# Define common physical and mathematical constants
${OBJDIR}/MyConsts.o : ${SRCDIR}/MyConsts.f90 ${OBJDIR}/ErrorTrap.o Makefile

# Set error and warning procedures
${OBJDIR}/ErrorTrap.o : ${SRCDIR}/ErrorTrap.f90 Makefile


