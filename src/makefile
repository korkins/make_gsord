#Comment: $make   $make install   $make clean
#https://www.gnu.org/software/make/manual/html_node/Suffix-Rules.html
.SUFFIXES: # erase current suffix list
.SUFFIXES: .f90 .o

F90SRC = AAA_TESTS_SORD_IP.f90 CSCATANG.f90 FORD_IP.f90 FORDM_IP.f90 FRESNELR0_IP.f90 FRESNELR_IP.f90 GAUSZW.f90 GROUPTAU.f90 NODES.f90 QX.f90\
         QX0.f90 RAYLEIGH12.f90 ROTATOR.f90 ROTATOR2.f90 ROUGHSRF.f90 RPV.f90 RTLS.f90 RTX.f90 SORD_IP.f90 SORDM_IP.f90\
	 SPLITTAU.f90 SSCAT1_IP.f90 SSCAT1M_IP.f90 SUMKA1B1.f90 SURFACE0_IP.f90 SURFACE0M_IP.f90 SURFACE_IP.f90 SURFACEM_IP.f90 SURFBNDM_IP.f90 SURFTOA_IP.f90\
	 Z0_IP.f90 ZAERM0_IP.f90 ZAERM_IP.f90 ZM0_IP.f90 ZM_IP.f90 ZRAYM0_IP.f90 ZRAYM_IP.f90

OFILES = ${F90SRC:.f90=.o}

FC=gfortran

SRCDIR = src

# HOME is system defined
INSTDIR = ${HOME}/mybins

COMPILEFLAGS = -O3

PROGRAM = sord

gsord.out : $(OFILES)
	${FC} $(OFILES) -o ${PROGRAM}
	
.f90.o:
	${FC} -c $(COMPILEFLAGS) $<

install:
	cp ${PROGRAM} ${INSTDIR}/${PROGRAM}
	
clean:
	rm *.o ${PROGRAM}
	
