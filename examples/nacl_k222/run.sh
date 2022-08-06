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
# 1:    0.0    0.0    0.0
# 2:    0.0    0.0    0.5
# 3:    0.0    0.5    0.0
# 4:    0.0    0.5    0.5
# 5:    0.5    0.0    0.0
# 6:    0.5    0.0    0.5
# 7:    0.5    0.5    0.0
# 8:    0.5    0.5    0.5

# $exe cell.txt basis_id.txt aims K 0 0 1/2 coulomb_matr_kpt_2.mtx

# $exe cell.txt basis_id.txt aims K 0 1/2 0 coulomb_matr_kpt_3.mtx

# $exe cell.txt basis_id.txt aims K 0 1/2 1/2 coulomb_matr_kpt_4.mtx

# $exe cell.txt basis_id.txt aims K 1/2 0 0 coulomb_matr_kpt_5.mtx

# $exe cell.txt basis_id.txt aims K 1/2 0 1/2 coulomb_matr_kpt_6.mtx

$exe cell.txt basis_id.txt aims K 1/2 1/2 0 coulomb_matr_kpt_7.mtx
# 0: 0.500000 0.500000 0.000000 -> 7
# 2: 0.000000 0.500000 0.500000 -> 4
# 4: 0.500000 0.000000 0.500000 -> 6

# $exe cell.txt basis_id.txt aims K 1/2 1/2 1/2 coulomb_matr_kpt_8.mtx
# 0: 0.500000 0.500000 0.500000 -> 8
# 2: 0.500000 0.000000 0.000000 -> 5
# 5: 0.000000 0.000000 0.500000 -> 2
# 6: 0.000000 0.500000 0.000000 -> 3
