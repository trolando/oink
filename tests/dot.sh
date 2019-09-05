#!/bin/bash
mkdir -p dot
for f in vb*
do
    if [ ! -e dot/$f.png ]; then
        ../build/dotty $f | dot -Tpng > dot/$f.png
    fi
#    if [ ! -e dot/$f-eo.png ]; then
#        ../build/nudge $f --evenodd | ../build/dotty | dot -Tpng > dot/$f.png
#    fi
done
