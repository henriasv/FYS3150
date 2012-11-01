#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=icc
CXX=icc
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/558680977/lib.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/Gauss_Laguerre.o \
	${OBJECTDIR}/Importance_Sampling_MC.o \
	${OBJECTDIR}/Bruteforce_MC.o \
	${OBJECTDIR}/Gauss_Legendre.o \
	${OBJECTDIR}/Helium_Solver.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-larmadillo -llapack -lblas -openmp -Wall
CXXFLAGS=-larmadillo -llapack -lblas -openmp -Wall

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/project3

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/project3: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/project3 ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/_ext/558680977/lib.o: ../lib/cppLibrary/lib.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/558680977
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/_ext/558680977/lib.o ../lib/cppLibrary/lib.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/Gauss_Laguerre.o: Gauss_Laguerre.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Gauss_Laguerre.o Gauss_Laguerre.cpp

${OBJECTDIR}/Importance_Sampling_MC.o: Importance_Sampling_MC.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Importance_Sampling_MC.o Importance_Sampling_MC.cpp

${OBJECTDIR}/Bruteforce_MC.o: Bruteforce_MC.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Bruteforce_MC.o Bruteforce_MC.cpp

${OBJECTDIR}/Gauss_Legendre.o: Gauss_Legendre.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Gauss_Legendre.o Gauss_Legendre.cpp

${OBJECTDIR}/Helium_Solver.o: Helium_Solver.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} $@.d
	$(COMPILE.cc) -g -MMD -MP -MF $@.d -o ${OBJECTDIR}/Helium_Solver.o Helium_Solver.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/project3

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
