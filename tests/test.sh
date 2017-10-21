#!/bin/bash
solvers=( "" "--pp" "--ppp" "--rr" "--dp" "--rrdp" "--psi -w -1" "--psi -w 0" "--zlk -w -1" "--zlk -w 0" "--qpt" "--spm" "--mspm" )
for SOLVER in "${solvers[@]}"
do
    echo "Testing solver $SOLVER"
    for f in vb*
    do
        # echo "Testing file $f"
        if [ "$SOLVER" = "--psi -w -1" ] || [ "$SOLVER" = "--psi -w 0" ]; then
            ../build/oink $f $SOLVER -v --no-loops > /dev/null
        else
            ../build/oink $f $SOLVER -v --no-loops --no-wcwc > /dev/null
        fi
        if [ $? -ne 0 ]; then
            echo "ERROR with solver $SOLVER and test file $f"
        fi
    done
done
