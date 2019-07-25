CC=mpic++
CFLAGS=-pedantic -Wall -Werror -Wextra -O3 -g -std=c++11 -Winline -fopenmp#-DNSERIALIZE
#INCLUDE=-I"$(BOOST_ROOT)"

source = $(wildcard *.cpp)
obj = $(source:.cpp=.o)
exe = $(basename $(source))
run_exe := $(run_$(basename $(source)))

all: clean $(exe)
	
$(obj): %.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	
$(exe): %: %.o
	$(CC) $(CFLAGS) -o $@ $<

run:
	mpiexec -n 16 ./advection 64 100 1.0 1.0e6 5.0e-7 2.85e-7 1 mpi_blocking

clean:
	rm -f $(obj) $(exe)