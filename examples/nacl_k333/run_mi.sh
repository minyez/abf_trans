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

# 0.000000000     0.000000000     0.000000000 11
# 0.500000000     0.500000000     0.500000000 17
# the basis set corresponding to a modified tier1 aims setup, i.e.
# Na
#   valence      3  s   1.
#   valence      2  p   6.
#   # hydro 2 p 1.2
#   # hydro 3 s 1.8
#    hydro 3 d 3.8
# Cl
#   valence      3  s   2.
#   valence      3  p   5.
#    ionic 3 d auto
#   # hydro 2 p 1.9
#    hydro 4 f 7.4
#   # ionic 3 s auto

# Note the following is the auxiliary basis, not the original basis
# ABFs on Na = 76, Cl = 197
cat > basis_id.txt << EOF
11 0 2
11 1 1
11 2 2
11 0 1
11 1 3
11 3 1
11 1 1
11 2 2
11 0 1
11 2 1
11 0 1
11 4 1
11 3 1
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
17 2 2
17 3 1
17 4 1
17 2 1
17 3 1
17 0 1
17 1 2
17 2 1
17 3 3
17 4 2
17 5 1
17 0 1
17 1 1
17 2 2
17 3 1
17 0 1
17 6 1
17 5 1
17 4 1
17 3 1
17 2 1
17 1 1
17 0 1
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

$exe cell.txt basis_id.txt matrix_inputs.txt

# $exe cell.txt basis_id.txt matrix_inputs_all.txt
