CC          =g++ -g
CFLAGS     = -O3 -fopenmp -Wall -I/opt/seas/pkg/gcc/fftw3/openmp/double/3.3.7/include
LIBS      = -lfftw3_omp -lfftw3 -lpthread -L/opt/seas/pkg/gcc/fftw3/openmp/double/3.3.7/lib


## CC 	   = mpic++ -std=c++11
#CC 	   = g++ -g 
## CC 	   = ~/usr/local/bin/mpic++ -std=c++11
##FFTW_LOC = /opt/seas/pkg/gcc/fftw3/mpi/double/3.3.7
#FFTW_LOC = ${HOME}/local/fftw3
## FFTWHOME = ${HOME}/local/fftw_mpi
## FFTW_LOC = include2/fftw_mpi
#EIGEN_LOC = include2/eigen
#EIGEN_PRINT_LOC = ${HOME}/Install/gdb
#INCLUDE = -I${FFTWHOME}/include -I/opt/seas/pkg/gcc/fftw3/openmp/double/3.3.7/include
#CFLAGS     =  -I${FFTWHOME}/include -I${EIGEN_LOC} -O3 -Wno-unused-result -Wno-write-strings  -fopenmp
## CFLAGS     = -I${FFTW_LOC}/include -I${FFTWHOME}/include -I${EIGEN_LOC} -O3 -Wno-unused-result -Wno-write-strings -fopenmp
#LIBS      = -lm  -lfftw3 -O3  -L${FFTWHOME}/lib   -L/opt/seas/pkg/gcc/fftw3/openmp/double/3.3.7/lib -lpthread -lfftw3_omp 
## LIBS      = -lm -lfftw3_mpi -lfftw3 -O3  -L${FFTWHOME}/lib -lfftw3_omp -lpthread  -L/opt/seas/pkg/gcc/fftw3/openmp/double/3.3.7/lib
#DLIBS      = -lm -lfftw3_mpi -lfftw3  -L${FFTWHOME}/lib
#DFLAGS     = -I${FFTW_LOC}/include  -I${FFTWHOME}/include -I${EIGEN_PRINT_LOC} -I${EIGEN_LOC} -Wno-unused-result -Wno-write-strings

#CC 	   = icpc 
#FFTW_LOC = ${HOME}/Install/fftw3a-threaded
#CFLAGS     = -O3 -openmp -I${FFTW_LOC}/include
#LIBS      = -openmp -O3 -lfftw3_threads -lfftw3 -lpthread -L${FFTW_LOC}/lib 

#############################################################################
# nothing should be changed below here

SRCS = template.cpp read_dump_traj.cpp grid_utils.cpp fftw_wrappers.cpp io_utils.cpp \
       die.cpp array_utils.cpp integ_utils.cpp random.cpp
			
			 
       
       
			 
DEBUGFOLDER := objects/debug
MAINOBJFOLDER := objects/main
DIM     = $(shell grep -e "^\#define Dim" globals.h | awk '{print $$3}')
IS2D = $(shell if [ $(grep -e "^\#define Dim" globals.h | awk '{print $$3}') == 2 ] echo true )
IS3D = $(shell if [[ $(DIM) == 3 ]] echo true )
# EXE := $(shell git describe --always)-$(DIM)d.out
TAGS := .tags

OBJS = ${SRCS:.cpp=.o}
DOBJS = $(addprefix $(DEBUGFOLDER)/, ${OBJS} )
OBJS_MAKE = $(addprefix $(MAINOBJFOLDER)/, ${OBJS} )


$(DEBUGFOLDER)/%.o: %.cpp
	@mkdir -p $(dir $@) 
	${CC} -g ${DFLAGS} -c $< -o $@  

$(MAINOBJFOLDER)/%.o: %.cpp
	@mkdir -p $(dir $@) 
	${CC} ${CFLAGS} ${DFLAGS} -c $< -o $@  

a_dmft-lc: ${OBJS_MAKE}
	$(CC) ${CFLAGS} ${DFLAGS} -o $@ ${OBJS_MAKE} $(LIBS)
ifeq ($(DIM), 2)
	@cp a_dmft-lc a_dmft-lc2d
endif 
ifeq ($(DIM), 3)
	@cp a_dmft-lc a_dmft-lc3d
endif 
	# ctags -R -f $(TAGS) .

# a_dmft-lc2d: k


debug: ${DOBJS}
	$(CC) -g  ${DFLAGS} -o $@ ${DOBJS} $(LIBS)

commit: a_dmft-lc
	cp a_dmft-lc $(EXE)

clean_debug: 
	rm -f $(DEBUGFOLDER)/*.o
	rm -f debug

debug_clean: 
	rm -f $(DEBUGFOLDER)/*.o
	rm -f debug
clean:
	rm -f $(MAINOBJFOLDER)/*.o
	rm -f a_dmft-lc
	rm -f *~

total_clean:
	rm -f a_dmft-lc* *~ debug *.dat LOG *.traj
	rm -fr objects

