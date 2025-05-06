ifdef TEC360HOME
CFLAGS		 = -I${TEC360HOME}/include/ -DTECIO=1 -g -O0
LIBS		 = ${TEC360HOME}/lib/tecio64.a -lstdc++
else
CFLAGS		 = -g -O0
LIBS             =
endif
LDFLAGS		 = 
FFLAGS		 =
CPPFLAGS 	 =  
FPPFLAGS         =
LOCDIR		 = 
MANSEC           = SNES

LIBFLAG          =

SOURCEC = compgeom.c ibm_iso.c ibm_io.c init.c\
          main.c metrics.c rhs.c wallfunction.c\
	  fsi.c solvers.c variables.c les.c  fish.c

OBJSC =  compgeom.o ibm_iso.o ibm_io.o init.o\
         main.o metrics.o rhs.o wallfunction.o\
         fsi.o solvers.o variables.o les.o fish.o

LIBBASE = libpetscmat


include ${PETSC_DIR}/lib/petsc/conf/variables

include ${PETSC_DIR}/lib/petsc/conf/rules

#include /sw/hprc/sw/petsc/3.6.2-intel-2017A-MPI-Hypr-debug/lib/petsc/conf/variables
#include /sw/hprc/sw/petsc/3.6.2-intel-2017A-MPI-Hypr-debug/lib/petsc/conf/rules



testt: ${OBJSC}
	-$(CLINKER) -o testt ${OBJSC} -O0 ${PETSC_LIB}

	rm main.o

data: variables.o compgeom.o data_ibm.o ibm_io.o wallfunction.o data.o
	-${CLINKER} -o data variables.o compgeom.o data_ibm.o ibm_io.o wallfunction.o data.o ${PETSC_LIB} ${LIBS}

itfcsearch: itfcsearch.o variables.o compgeom.o
	-${CLINKER} -o itfcsearch itfcsearch.o  variables.o compgeom.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB}

data_vtk: data_surface.o 
	-${CLINKER} -o data_vtk data_surface.o  ${PETSC_SNES_LIB} ${PETSC_TS_LIB}

datalis: data_file2lis.o
	-${CLINKER} -o datalis data_file2lis.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}

datafile: data_list2file.o
	-${CLINKER} -o datafile data_list2file.o ${PETSC_SNES_LIB} ${PETSC_TS_LIB} ${LIBS}
gridgen: gridgen.o
	-${CLINKER} -o gridgen gridgen.o ${PETSC_LIB} ${LIBS}
cleanobj:
	rm -f *.o

include ${PETSC_DIR}/lib/petsc/conf/test
#include /sw/hprc/sw/petsc/3.6.2-intel-2017A-MPI-Hypr-debug/lib/petsc/conf/test

