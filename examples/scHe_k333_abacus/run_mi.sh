#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

cat > cell.txt << EOF
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
1
0.000000000     0.000000000     0.000000000 2
EOF

# basis on He = 14
cat > basis_id.txt << EOF
2 0 2
2 1 1
EOF

$exe cell.txt basis_id.txt matrix_inputs.txt
