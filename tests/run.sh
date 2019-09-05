#!/bin/bash
SOLVER=$@
for f in ../tests/vb*
do
    LOGFILE=`mktemp` || exit 1
    ../build/oink $f $SOLVER --no-loops --no-wcwc -t -v > $LOGFILE
    if [ $? -ne 0 ]; then
        echo "ERROR: $f"
    else
        DUPLICATES=`grep -o "duplicate" $LOGFILE | wc -l`
        SUBS=`grep -o "subtangle" $LOGFILE | wc -l`
        echo "GOOD: $f ($DUPLICATES duplicates, $SUBS subtangles)"
    fi
    rm -f $LOGFILE
done
