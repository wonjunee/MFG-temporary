#!/bin/sh

# CC := g++
CC := g++ # typical homebrew compiler name 

# Check for environment definitions of compiler 
# e.g., on CC = g++-6 on OSX
ifdef PHYSICELL_CPP 
	CC := $(PHYSICELL_CPP)
endif

 
CFLAGS     :=  -O3 -std=c++17 -I./lib
# FFTW library
FFTWFLAG   := -I/usr/local/include  -L/usr/local/lib -lfftw3
# opencv library: for image
OPENCVFLAG := `pkg-config --cflags --libs opencv4`
# for parallelization. Use 0 if you don't want. Use the multiple of 2 if possible.
OPENMPFLAG := -DASYNC=10
 
COMPILE_COMMAND := $(CC) $(CFLAGS) $(FFTWFLAG) $(OPENCVFLAG) $(OPENMPFLAG)
 
OUTPUT := MFG
 
all: ./lib/main.cpp
	$(COMPILE_COMMAND) -o $(OUTPUT) ./lib/main.cpp

clean:
	rm -f *.o $(OUTPUT).*
