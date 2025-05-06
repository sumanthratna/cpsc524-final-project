#!/usr/bin/env bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 <prog> <graph> <threshold>"
    exit 1
fi

GRAPH=$1
THRESH=$2

# Get number of vertices based on graph name
case "$GRAPH" in
    *"amazon0601"*) V=403394 ;;
    *"gplus_combined"*) V=107614 ;;
    *"soc-pokec-relationships"*) V=1632804 ;;
    *"twitter_combined"*) V=81306 ;;
    *"web-BerkStan"*) V=685230 ;;
    *"web-Google"*) V=916428 ;;
    *"web-NotreDame"*) V=325729 ;;
    *"web-Stanford"*) V=281904 ;;
    *"wiki-Talk"*) V=2394385 ;;
    *"soc-Epinions1"*) V=76000 ;;
    *) echo "Unknown graph: $GRAPH"; exit 1 ;;
esac

C_PROGS=("serial_jacobi" "jacobi_blocking" "guess_seidel" "guess_seidel_cache" "gauss_seidel_blocking")
CXX_PROGS=("delta_push" "parlay_jacobi" "random_walk")

echo "Running all programs 10 times..."

for prog in "${C_PROGS[@]}" "${CXX_PROGS[@]}"; do
    echo "Running ./$prog $GRAPH $V $THRESH 0.85 8 ten times..."
    for i in {1..10}; do
        echo "Run #$i for $prog"
        ./$prog $GRAPH $V $THRESH 0.85 8
    done
done
