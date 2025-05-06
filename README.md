# cpsc524-final-project

Sumanth Ratna, Micah Gold, Buwei Chen

## Implementations

- [`serial_jacobi`](./serial_jacobi.c)
- [`jacobi_blocking`](./jacobi_blocking.c)
- [`guess_seidel`](./guess_seidel.c)
- [`guess_seidel_cache`](./guess_seidel_cache.c)
- [`gauss_seidel_blocking`](./gauss_seidel_blocking.c)
- [`delta_push`](./delta_push.cpp)
- [`parlay_jacobi`](./parlay_jacobi.cpp)
- [`random_walk`](./random_walk.cpp)

## Instructions to Run

Potentially need to change g++ location (line 8-9) in makefile.

1. To compile: `make`

2. To run one implementation: `./run_benchmark.sh delta_push datasets/web-Google.txt 0.000001`

You can see a list of datasets [here](./datasets/). The third parameter is the convergence threshold. 
