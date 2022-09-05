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

# ABFs on He = 14
cat > basis_id.txt << EOF
2 0 3
2 1 2
2 2 1
EOF

rm -f abf_trans_out*

$exe cell.txt basis_id.txt matrix_inputs.txt

tested=$(tail -4 log.txt | awk '/Tested transformations/ {print $3}')
passed=$(tail -4 log.txt | awk '/Passed transformations/ {print $3}')
(( tested != passed )) && exit 1
exit 0
