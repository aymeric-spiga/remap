#CFLAGS = -ansi -pedantic -w -g -O2 -fPIC #-DDEBUG $(EXTRA_FLAGS)
#CFLAGS = -ansi -g -w -O2 -fPIC
#CFLAGS = -ansi -w -O0 -fPIC
CFLAGS = -ansi -w -O3 -fPIC
#CXX = OMPI_CXX=/opt/intel/bin/icpc mpic++ -DOMPI_SKIP_MPICXX
CXX = /u/aslmd/miniconda/bin/mpicxx -DOMPI_SKIP_MPICXX
CC = /u/aslmd/miniconda/bin/mpicc

vpath %.cpp src
SRC = circle.cpp cputime.cpp grid.cpp inside.cpp intersect.cpp \
	mapper.cpp meshutil.cpp mpi_routing.cpp mpi_cascade.cpp \
	node.cpp parallel_tree.cpp polyg.cpp \
	timer.cpp tree.cpp triple.cpp \
	libmapper.cpp

OBJ = $(patsubst %.cpp,obj/%.o,$(SRC))

lib/libmapper.so: $(OBJ)
	mkdir -p lib
	$(CXX) -shared -o $@ $^ -lc

test: $(OBJ) obj/test-main.o
	$(CXX) -o $@ $^ -lnetcdf

try: test
	cd run && mpirun -n 32 ./test

obj/%.o: %.cpp
	@mkdir -p obj
	$(CXX) -c $< -o $@ $(CFLAGS)

clean:
	rm -rf obj/*
