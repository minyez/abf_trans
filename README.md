---
title: Demo of calculation of transformation matrix between equivalent k/R points
---

This project briefly implements the transformation matrix to obtain operator representation at one k/R point
to other equivalent points within the atom-centered basis functions framework.

## Compilation

1. Ensure a valid `SPGLIB_HOME` environment variable as install path of [Spglib](https://spglib.github.io).
2. Ensure a valid `LAPACK_HOME` environment variable and correct librarires in `LAPACK_LIBRARIES` in CMakeLists.txt.
   For example, if you use MKL,
    ```
    # in bashrc or envrc
    export LAPACK_HOME="$MKLROOT/lib/intel64"
    # in CMakeLists.txt
    set(LAPACK_LIBRARIES "-L$ENV{LAPACK_PATH} -lmkl_intel_lp64 -lmkl_sequential -lmkl_core")
    ```
3. (optional) Customize `CMAKE_CXX_COMPILER` and `CMAKE_CXX_FLAGS` in CMakeLists.txt to suit your need.
4. Run the general cmake process

    ```bash
    mkdir build && cd && build
    cmake .. && make -j4
    ```

This will leads to executable `abf_trans.exe`, and test executables in `build/tests`.
Run `make test` to do some basic unit tests.

## Usage

### Inputs

Two files are required as inputs to use this project

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
   type_a l_rf
   ```
  where `type_a` is the type of atom and `l_rf` the angular momentum quantum number of the radial function.
  There can be lines that are identical to each other, but the radial functions are different by definition.

To run the program
```shell
abf_trans.exe cell.txt basis_id.txt choice mode x1 x2 x3 mtxfile
```
where `choice` is the choice of real spherical harmonics. It can be `orig` or `aims` for the current stage.
`mode` should be either `R`/`K`,
and `xN` are the components along direct/reciprocal lattice vectors.
`mtxfile` is the name of matrix data file in sparse matrix-market format for the matrix representation at the specified k/R-point.

Another way to run the program is

```shell
abf_trans.exe cell.txt basis_id.txt matrices.txt
```
where each line in `matrices.txt` contain the `choice mode x1 x2 x3 mtxfile` information.
   
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
```
