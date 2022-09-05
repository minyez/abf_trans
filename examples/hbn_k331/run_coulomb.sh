#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

cat > cell.txt << EOF
 2.510000    0.000000    0.000000
-1.255000    2.173724    0.000000
 0.000000    0.000000    6.690000
4
0.0000     0.000     0.250000000  7
0.0000     0.000    -0.250000000  7
1./3.0     -1./3.    0.250000000  5
-1./3.0     1./3.   -0.250000000  5
EOF

cat > basis_id.txt << EOF
7 0 2
7 1 2
7 0 1
7 2 1
7 1 1
7 0 1
5 0 2
5 1 1
5 0 1
5 1 1
5 2 1
5 1 1
5 0 1
EOF

$exe cell.txt basis_id.txt matrix_inputs.txt
tested=$(tail -4 log.txt | awk '/Tested transformations/ {print $3}')
passed=$(tail -4 log.txt | awk '/Passed transformations/ {print $3}')
(( tested != passed )) && exit 1
exit 0
