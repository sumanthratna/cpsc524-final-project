# C Compiler + Flags
CC = gcc
CFLAGS = -O4
CLIBS = -lpthread -lm

# C++ Compiler + Flags
CXX = /opt/homebrew/opt/gcc/bin/g++-14 # for buwei
# CXX = g++ # for micah
CXXFLAGS = -std=c++17 -O3 -fopenmp -pthread
PARLAY_INC = -I parlaylib/include

# Executables (C and C++)
C_PROGS = serial_jacobi jacobi_blocking guess_seidel guess_seidel_cache gauss_seidel_blocking
CXX_PROGS = delta_push parlay_jacobi random_walk gmres

# All executables
PROGS = $(C_PROGS) $(CXX_PROGS)

all: $(PROGS)

# C rule
$(C_PROGS): %: %.c
	$(CC) $(CFLAGS) $< -o $@ $(CLIBS)

# C++ rule (uses Parlay include for all C++ sources)
$(CXX_PROGS): %: %.cpp
	$(CXX) $(CXXFLAGS) $(PARLAY_INC) $< -o $@

# Benchmark single
benchmark:
		EXEC=$(word 2,$(MAKECMDGOALS)); \
		GRAPH=$(word 3,$(MAKECMDGOALS)); \
		THRESH=$(word 4,$(MAKECMDGOALS)); \
		echo "Running ./$$EXEC $$GRAPH 916428 $$THRESH 0.85 8"; \
		./$$EXEC $$GRAPH 916428 $$THRESH 0.85 8; \
	exit 0

# Benchmark all
benchmark-all:
	GRAPH=$(word 2,$(MAKECMDGOALS)); \
	THRESH=$(word 3,$(MAKECMDGOALS)); \
	for prog in $(PROGS); do \
		echo "=== Running $$prog with $$GRAPH $$THRESH ==="; \
		./$$prog $$GRAPH 916428 $$THRESH 0.85 8 || echo "$$prog failed"; \
		echo ""; \
	done; exit 0

clean:
	rm -f $(PROGS) *.o *.out *.exe *.bin *.txt

# Hack to ignore extra args
%::
	@:
