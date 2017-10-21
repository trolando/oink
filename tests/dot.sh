#!/bin/bash
mkdir dot
for f in vb*
do
    ../build/dotty $f | dot -Tpng > dot/$f.png
done
