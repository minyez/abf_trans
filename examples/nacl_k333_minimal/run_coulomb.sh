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

# Na
#   valence      3  s   1.
#   valence      2  p   6.
# Cl
#   valence      3  s   2.
#   valence      3  p   5.

# Note the following is the auxiliary basis, not the original basis
# ABFs on Na = 76, Cl = 197
cat > basis_id.txt << EOF
11 0 2
11 1 1
11 0 1
11 1 2
11 0 1
11 2 1
11 1 1
11 0 1
17 0 2
17 1 1
17 0 1
17 1 2
17 2 1
17 1 1
17 0 2
17 1 1
17 2 1
17 0 1
17 1 1
17 2 1
17 1 1
17 0 1
EOF

rm -f abf_trans_out*

$exe cell.txt basis_id.txt matrix_inputs.txt

tested=$(tail -4 log.txt | awk '/Tested transformations/ {print $3}')
passed=$(tail -4 log.txt | awk '/Passed transformations/ {print $3}')
(( tested != passed )) && exit 1
exit 0
