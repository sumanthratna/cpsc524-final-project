#!/usr/bin/env bash

if [ $# -lt 3 ]; then
    echo "Usage: $0 <prog> <graph> <threshold>"
    exit 1
fi

EXEC=$1
GRAPH=$2
THRESH=$3

# Get number of vertices based on graph name
case "$GRAPH" in
    *"amazon0601"*) V=403394 ;;
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

echo "Running ./$EXEC $GRAPH $V $THRESH 0.85 8"
./$EXEC $GRAPH $V $THRESH 0.85 8
