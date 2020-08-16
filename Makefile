# Make file for building application

CC = g++
CFLAGS = -Wall -Wextra -std=c++11 -lpthread -pedantic -O3 -ldl -pthread -fopenmp

all: volume simd

volume: volume.o
	$(CC) volume.o $(CFLAGS) -o volume

volume.o: volume.cpp volume.hpp
	$(CC) -c volume.cpp $(CFLAGS)

simd: simd.o
	$(CC) simd.o $(CFLAGS) -o simd

simd.o: simd.cpp simd.hpp
	$(CC) -c simd.cpp $(CFLAGS)

clean:
	rm -f *.o *.d *.gch volume simd
