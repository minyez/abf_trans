#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

cat >> cell.txt << EOF
 2.510000    0.000000    0.000000
-1.255000    2.173724    0.000000
 0.000000    0.000000    6.690000
4
0.000000000     0.000000000     0.250000000  7
0.000000000     0.000000000     0.750000000  7
0.333333333     0.666666667     0.250000000  5
0.666666667     0.333333333     0.750000000  5
EOF

cat >> basis_id.txt << EOF
5 1 0
5 2 0
5 3 1
7 4 0
7 5 0
7 6 1
EOF

$exe cell.txt basis_id.txt K 1/2 1/2 1/2
