#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

cat > cell.txt << EOF
0.00000000    2.82865000    2.82865000
2.82865000    0.00000000    2.82865000
2.82865000    2.82865000    0.00000000
2
0.000000000     0.000000000     0.000000000 11
0.500000000     0.500000000     0.500000000 17
EOF

# the basis set corresponding to a tier1 aims setup, i.e.
# Na
#   valence      3  s   1.
#   valence      2  p   6.
#    hydro 2 p 1.2
#    hydro 3 s 1.8
#    hydro 3 d 3.8
# Cl
#   valence      3  s   2.
#   valence      3  p   5.
#    ionic 3 d auto
#    hydro 2 p 1.9
#    hydro 4 f 7.4
#    ionic 3 s auto

cat > basis_id.txt << EOF
11 0
11 0
11 0
11 0
11 1
11 1
11 2
17 0
17 0
17 0
17 0
17 1
17 1
17 1
17 2
17 3
EOF

$exe cell.txt basis_id.txt K 1/2 1/2 1/2
