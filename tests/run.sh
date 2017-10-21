#!/bin/bash
SOLVER=$1
for f in vb*
do
    ../build/oink $f $SOLVER --no-loops --no-wcwc -v > /dev/null
    if [ $? -ne 0 ]; then
        echo "ERROR: $f"
    else
        echo "GOOD: $f"
    fi
done
