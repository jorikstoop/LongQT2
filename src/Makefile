# General Code Structure (GCS) for Makefile
# 
# Under Intellectual Property Protection.
# 
# Author      	: Shugang Zhang
# Modified by	: Jorik Stoop 
# Last update	: Nov-2026

all: SingleCell OneD Initialization

common = TP06.cc  Cell.cc

CC = g++
CFLAGS = -std=c++17 -O3 -w #-g:warning output to screen   -w:warning ignored
OMPFLAGS = -fopenmp

SingleCell: $(common) SingleCell.cc
	$(CC) $(CFLAGS) -o model_single_cell $(common) SingleCell.cc

Initialization: $(common) Initialization.cc
	$(CC) $(CFLAGS) -o model_initialization $(common) Initialization.cc

OneD: $(common) OneD.cc
	$(CC) $(CFLAGS) $(OMPFLAGS) -o model_oned $(common) OneD.cc

clean:
	del model_*

# Note: use rm model_* for Unix but del model_* for windows
