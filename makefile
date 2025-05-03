CC = gcc
FLG = -O4
LDLIBS = -lpthread -lm

PROGS = serial_jacobi jacobi_blocking guess_seidel guess_seidel_batched guess_seidel_cache gauss_seidel_blocking

all: $(PROGS)

%: %.c
	$(CC) $(FLG) $< -o $@ $(LDLIBS)

# Positional benchmark: make benchmark <prog> <graph> <threshold>
benchmark:
	@if [ "$(word 1,$(MAKECMDGOALS))" = "benchmark" ]; then \
		echo "Usage: make benchmark <prog> <graph> <threshold>"; \
	else \
		EXEC=$(word 2,$(MAKECMDGOALS)); \
		GRAPH=$(word 3,$(MAKECMDGOALS)); \
		THRESH=$(word 4,$(MAKECMDGOALS)); \
		echo "Running ./$$EXEC $$GRAPH 916428 $$THRESH 0.85 8"; \
		./$$EXEC $$GRAPH 916428 $$THRESH 0.85 8; \
	fi; exit 0

# Positional benchmark-all: make benchmark-all <graph> <threshold>
benchmark-all:
	GRAPH=$(word 2,$(MAKECMDGOALS)); \
	THRESH=$(word 3,$(MAKECMDGOALS)); \
	for prog in $(PROGS); do \
		echo "=== Running $$prog with $$GRAPH $$THRESH ==="; \
		./$$prog $$GRAPH 916428 $$THRESH 0.85 8 || echo "$$prog failed"; \
		echo ""; \
	done; exit 0

clean:
	rm -f $(PROGS) *.o *.out *.exe *.bin

# Hack to prevent Make from treating args as targets
%::
	@:
