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
LDFLAGS  += -L/opt/homebrew/opt/libomp/lib -lomp
OMPFLAGS = -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include

SingleCell: $(common) SingleCell.cc
	$(CC) $(CFLAGS) $(OMPFLAGS) $(LDFLAGS) -o model_single_cell $(common) SingleCell.cc

Initialization: $(common) Initialization.cc
	$(CC) $(CFLAGS) $(OMPFLAGS) $(LDFLAGS) -o model_initialization $(common) Initialization.cc

OneD: $(common) OneD.cc
	$(CC) $(CFLAGS) $(OMPFLAGS) $(LDFLAGS) -o model_oned $(common) OneD.cc

clean:
	rm model_*

# Note: use rm model_* for Unix but del model_* for windows
