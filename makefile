# Compiler settings
CC = gcc
# CXX = /opt/homebrew/opt/gcc/bin/g++-14 # for buwei
CXX = g++
CFLAGS = -O4
CXXFLAGS = -std=c++17 -O3 -fopenmp -pthread
CLIBS = -lpthread -lm
PARLAY_INC = -I parlaylib/include

# Program lists
C_PROGS = serial_jacobi jacobi_blocking guess_seidel guess_seidel_cache gauss_seidel_blocking
CXX_PROGS = delta_push parlay_jacobi random_walk
PROGS = $(C_PROGS) $(CXX_PROGS)

# Default target
all: $(PROGS)

# Build rules
$(C_PROGS): %: %.c
	$(CC) $(CFLAGS) $< -o $@ $(CLIBS)

$(CXX_PROGS): %: %.cpp
	$(CXX) $(CXXFLAGS) $(PARLAY_INC) $< -o $@

# Clean target
clean:
	rm -f $(PROGS) *.o *.out *.exe *.bin

.PHONY: all clean
