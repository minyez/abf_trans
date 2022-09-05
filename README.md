# Symmetry transformation of operator representations on equivalent k/R points

This project briefly implements the transformation of operator representation
in atom-centered basis functions at one k/R point
to its counterparts at equivalent points determined by the symmetry of the system.

**N.B.**: The project is for demonstration and practice purpose.
At the current stage, only k-space symmetry transformation is implemented.

## Requirement

- C++ compiler with C++-11 support
- [Spglib](https://spglib.github.io)
- A linear-algebra package, e.g. [LAPACK](https://netlib.org/lapack/)

## Compilation

1. Ensure a valid `SPGLIB_HOME` environment variable as install path of Spglib.
2. Ensure a valid `LAPACK_HOME` environment variable and correct librarires in `LAPACK_LIBRARIES` in CMakeLists.txt.
   For example, if using MKL,
    ```
    # in bashrc or envrc
    export LAPACK_HOME="$MKLROOT/lib/intel64"
    # in CMakeLists.txt
    set(LAPACK_LIBRARIES "-L$ENV{LAPACK_PATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
    ```
3. (optional) Customize `CMAKE_CXX_COMPILER` and `CMAKE_CXX_FLAGS` in CMakeLists.txt to suit your need.
4. Run the general cmake procedure

    ```bash
    mkdir build && cd && build
    cmake .. && make -j4
    ```

This will lead to executable `abf_trans.exe` under `build` directory and test executables in `build/tests`.
You may run `make test` to do basic unit tests.

## Usage

The main program `abf_trans.exe` can be run by the following command

```shell
abf_trans.exe cell.txt basis_id.txt matrix_inputs
```

After start, the program will
1. read the cell structure and generate space-group symmetry operations.
2. read the basis and setup the indices including radial functions, angular and azimuth quantum numbers.
3. read matrix representations at different k/R points.
4. for each k/R point,
    1. find its equivalent point using the symmetry operations
    2. perform the unitary transformation of the matrix, and compare the transformed matrix to the one at the equivalent point.

Three input files are required as arguments to use this project

- `cell.txt` for parsing the information of the cell.
    ```
    a_11 a_12 a_13
    a_21 a_22 a_23
    a_31 a_32 a_33
    n_atoms
    x1_a1 x2_a1 x3_a1 type_a1
    x1_a2 x2_a2 x3_a2 type_a2
    ...
    x1_an x2_an x3_an type_an
    ```
    where `xN` are the fractional coordinate of atoms and `type` (integer) is the identifier of the atom type.
- `basis_id.txt` to identify the information of basis functions.
  Each line of the file represents a set of basis functions for one atom species with the same radial function
  The format is
   ```
   type_a l_rf n_l
   ```
  where `type_a` is the type of atom, `l_rf` the angular momentum quantum number of the radial functions, and `n_l` is the number of different radial functions in this channel.
  There can be lines that have identical `type_a` and `l_rf`, but with a different `n_l`, for example, `n_l1` and `n_l2`.
  In this case, `n_l1` plus `n_l2` radial functions in the `l_rf` channel is added to `type_a` specie.
- `matrix_inputs` file is used to specify the data and related information of matrix representations to transform.
  It has the following format
  ```
  choice mode ngrid1 ngrid2 ngrid3
  x1_1 x2_1 x3_1 matfile_1
  x1_2 x2_2 x3_2 matfile_2
  # ...
  x1_n x2_n x3_n matfile_n
  ```
  In the header, `choice` is the choice of conventions, e.g. real spherical harmonics, Bloch sum, etc.
  They are generally chosen by the code used to generate the matrices.
  At the current stage, only the convetions of FHI-aims (`aims`) and ABACUS (`abacus`) are implemented.
  `mode` should be either `R`/`K`. `ngrid1`, `ngrid2` and `ngrid3` are the numbers of k/R grids on the three directions.
  For each of the `n` lines below the header,  `x1_n`, `x2_n` and `x3_n` are the coordinates of k (in reciprocal lattice vectors)  or R (in direct lattice vectors) point.
  `matfile_n` is the name of data file, containing the matrix representation at the specified k/R-point.
  File format of sparse matrix-market format (`mtx`) or ELSI CSC file (`csc`) is supported.
  The code will check if the specified matrix files exist and whether the contained data are valid.

Examples are given in the `examples` directory.
Please check `matrix_inputs` files and shell scripts under the subdirectories therein.

## References

```bibtex
@article{BlancoMA97,
  title = {Evaluation of the Rotation Matrices in the Basis of Real Spherical Harmonics},
  author = {Blanco, Miguel A. and Fl\'orez, M. and Bermejo, M.},
  date = {1997-12-08},
  journaltitle = {J. Mol. Struct. THEOCHEM},
  volume = {419},
  number = {1},
  pages = {19--27},
  issn = {0166-1280},
  doi = {10.1016/S0166-1280(97)00185-1},
  url = {https://www.sciencedirect.com/science/article/pii/S0166128097001851},
}

@article{DovesiR86,
  title = {On the Role of Symmetry in the Ab Initio Hartree-Fock Linear-Combination-of-Atomic-Orbitals Treatment of Periodic Systems},
  author = {Dovesi, Roberto},
  date = {1986},
  journaltitle = {Int. J. Quantum Chem.},
  volume = {29},
  number = {6},
  pages = {1755--1774},
  issn = {1097-461X},
  doi = {10.1002/qua.560290608},
  url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/qua.560290608},
}

@article{ZhangZ22,
  title = {{{MagneticTB}}: {{A}} Package for Tight-Binding Model of Magnetic and Non-Magnetic Materials},
  author = {Zhang, Zeying and Yu, Zhi-Ming and Liu, Gui-Bin and Yao, Yugui},
  date = {2022-01},
  journaltitle = {Comput. Phys. Commun.},
  volume = {270},
  pages = {108153},
  issn = {00104655},
  doi = {10.1016/j.cpc.2021.108153},
  url = {https://linkinghub.elsevier.com/retrieve/pii/S0010465521002654},
}
```
