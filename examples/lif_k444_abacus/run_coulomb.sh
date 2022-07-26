#!/bin/bash

exe="../../build/abf_trans.exe"

[[ ! -x $exe ]] && echo "Executable not built in ../../build/" && exit 1

cat > cell.txt << EOF
0.00000000    2.01350000    2.01350000
2.01350000    0.00000000    2.01350000
2.01350000    2.01350000    0.00000000
2
0.000000000     0.000000000     0.000000000 3
0.500000000     0.500000000     0.500000000 9
EOF

# Li  Li_gga_9au_100Ry_4s1p.orb
# F   F_gga_7au_100Ry_2s2p1d.orb

cat > basis_id.txt << EOF
3 0 6
3 1 4
3 2 1
9 0 5
9 1 4
9 2 4
9 3 2
9 4 1
EOF

rm -f abf_trans_out*

# the coordinates of corresponding k indices are

#  1:  0.0   0.0   0.0 
#  2:  0.0   0.0   1/3 
#  3:  0.0   0.0   2/3 
#  4:  0.0   1/3   0.0 
#  5:  0.0   1/3   1/3 
#  6:  0.0   1/3   2/3 
#  7:  0.0   2/3   0.0 
#  8:  0.0   2/3   1/3 
#  9:  0.0   2/3   2/3 
# 10:  1/3   0.0   0.0 
# 11:  1/3   0.0   1/3 
# 12:  1/3   0.0   2/3 
# 13:  1/3   1/3   0.0 
# 14:  1/3   1/3   1/3 
# 15:  1/3   1/3   2/3 
# 16:  1/3   2/3   0.0 
# 17:  1/3   2/3   1/3 
# 18:  1/3   2/3   2/3 
# 19:  2/3   0.0   0.0 
# 20:  2/3   0.0   1/3 
# 21:  2/3   0.0   2/3 
# 22:  2/3   1/3   0.0 
# 23:  2/3   1/3   1/3 
# 24:  2/3   1/3   2/3 
# 25:  2/3   2/3   0.0 
# 26:  2/3   2/3   1/3 
# 27:  2/3   2/3   2/3

$exe cell.txt basis_id.txt ./matrix_inputs_coulomb.txt

tested=$(tail -4 log.txt | awk '/Tested transformations/ {print $3}')
passed=$(tail -4 log.txt | awk '/Passed transformations/ {print $3}')
(( tested != passed )) && exit 1
exit 0
