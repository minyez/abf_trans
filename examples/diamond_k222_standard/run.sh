#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

# the first atom should direct mapped to (0,0,1/2) in the standard cell by transformation matrix
# but there is still a non-zero origin shift. Not quite clear about the mechanism
cat > cell.txt << EOF
0.000000000     1.783560000     1.783560000
1.783560000     0.000000000     1.783560000
1.783560000     1.783560000     0.000000000
2
0.500000000     0.500000000    -0.500000000   6
0.750000000     0.750000000    -0.250000000   6
EOF

# the basis set corresponding to a tier1 aims setup, i.e.
#   valence      2  s   2.
#   valence      2  p   2.
#   hydro 2 p 1.7
#   hydro 3 d 6
#   hydro 2 s 4.9

cat > basis_id.txt << EOF
6 0
6 0
6 0
6 1
6 1
6 2
EOF


$exe cell.txt basis_id.txt K 1/2 1/2 1/2
