#!/bin/bash
mkdir -p dot
for f in vb*
do
    if [ ! -e dot/$f.png ]; then
        ../build/dotty $f | dot -Tpng > dot/$f.png
    fi
done
