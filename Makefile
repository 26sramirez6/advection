CC=mpic++
CFLAGS=-pedantic -Wall -Werror -Wextra -O3 -g -std=c++11 -Winline -fopenmp #-DNSERIALIZE

source = $(wildcard *.cpp)
obj = $(source:.cpp=.o)
exe = $(basename $(source))

all: clean $(exe)
	
$(obj): %.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@
	
$(exe): %: %.o
	$(CC) $(CFLAGS) -o $@ $<

benchmarks: strong_blocking strong_non_blocking weak_blocking weak_non_blocking hybrid 

serial:
	mpiexec -n 1 ./advection 64 100 1.0 1.0e6 5.0e-7 2.85e-7 1 serial
	python main.py 64 1 serial
	
blocking:
	mpiexec -n 16 ./advection 64 100 1.0 1.0e6 5.0e-7 2.85e-7 1 mpi_blocking
	python main.py 64 1 blocking

non_blocking:
	mpiexec -n 16 ./advection 64 100 1.0 1.0e6 5.0e-7 2.85e-7 1 mpi_non_blocking
	python main.py 64 1 non_blocking

strong_blocking:
	rm -f strong_blocking.txt
	mpiexec -n 1 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 4 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 9 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 16 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 25 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 36 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 49 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 64 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 81 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt
	mpiexec -n 100 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> strong_blocking.txt

strong_non_blocking:
	rm -f strong_non_blocking.txt
	mpiexec -n 1 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 4 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 9 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 16 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 25 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 36 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 49 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 64 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 81 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt
	mpiexec -n 100 ./advection 5040 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> strong_non_blocking.txt

weak_blocking:
	rm -f weak_blocking.txt
	mpiexec -n 1 ./advection 800 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 4 ./advection 1600 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 9 ./advection 2400 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 16 ./advection 3200 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 25 ./advection 4000 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 36 ./advection 4800 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 49 ./advection 5600 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 64 ./advection 6400 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 81 ./advection 7200 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt
	mpiexec -n 100 ./advection 8000 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_blocking >> weak_blocking.txt

weak_non_blocking:
	rm -f weak_non_blocking.txt
	mpiexec -n 1 ./advection 800 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 4 ./advection 1600 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 9 ./advection 2400 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 16 ./advection 3200 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 25 ./advection 4000 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 36 ./advection 4800 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 49 ./advection 5600 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 64 ./advection 6400 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 81 ./advection 7200 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt
	mpiexec -n 100 ./advection 8000 300 1.0 1.0e3 5.0e-7 2.85e-7 1 mpi_non_blocking >> weak_non_blocking.txt

hybrid:
	rm -f hybrid.txt
	mpiexec -n 1 ./advection 10000 200 1.0 1.0e3 5.0e-7 2.85e-7 16 hybrid >> hybrid.txt
	mpiexec -n 4 ./advection 10000 200 1.0 1.0e3 5.0e-7 2.85e-7 4 hybrid >> hybrid.txt
	mpiexec -n 16 ./advection 10000 200 1.0 1.0e3 5.0e-7 2.85e-7 1 hybrid >> hybrid.txt

clean:
	rm -f $(obj) $(exe)