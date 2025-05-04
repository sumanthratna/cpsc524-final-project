#!/usr/bin/env bash

if [ $# -lt 3 ]; then
    echo "Usage: $0 <prog> <graph> <threshold>"
    exit 1
fi

EXEC=$1
GRAPH=$2
THRESH=$3

echo "Running ./$EXEC $GRAPH 916428 $THRESH 0.85 8"
./$EXEC $GRAPH 916428 $THRESH 0.85 8
