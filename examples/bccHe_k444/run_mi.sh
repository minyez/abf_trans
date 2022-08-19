#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

cat > cell.txt << EOF
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
2
0.000000000     0.000000000     0.000000000 2
0.500000000     0.500000000     0.500000000 2
EOF

# He
#   valence      1  s   2.
#   hydro        2  p   3.5

# ABFs on He = 13
cat > basis_id.txt << EOF
2 1 1
2 0 1
2 2 1
2 1 1
2 0 1
EOF

rm -f abf_trans_out*

$exe cell.txt basis_id.txt matrix_inputs.txt
