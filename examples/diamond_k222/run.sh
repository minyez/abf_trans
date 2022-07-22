#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

cat >> cell.txt << EOF
0.000000000     1.783560000     1.783560000
1.783560000     0.000000000     1.783560000
1.783560000     1.783560000     0.000000000
2
0.000000000     0.000000000     0.000000000   6
0.250000000     0.250000000     0.250000000   6
EOF

# the basis set corresponding to a tier1 aims setup, i.e.
#   valence      2  s   2.
#   valence      2  p   2.
#   hydro 2 p 1.7
#   hydro 3 d 6
#   hydro 2 s 4.9

cat >> basis_id.txt << EOF
6 0
6 0
6 0
6 1
6 1
6 2
EOF


$exe cell.txt basis_id.txt K 1/2 1/2 1/2
