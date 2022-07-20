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

This will leads to executable `abf_rotation.exe`, and test executables in `build/tests`.
Run `make test` to do some basic unit tests.

## Usage

### Inputs

Three files are required as inputs to use this project

- `lattic.txt` for parsing the information of the lattice.
    ```
    a_11 a_12 a_13
    a_21 a_22 a_23
    a_31 a_32 a_33
    n_atoms
    x_a1 y_a1 z_a1 type_a1
    x_a2 y_a2 z_a2 type_a2
    ...
    x_an y_an z_an type_an
    ```
- `basis_id.txt` to identify the information of basis functions.
  Each line of the file represents a set of basis functions for one atom specie, with the following format
   ```
   type_a t n l
   ```
  where `type_a` is the type of atom, `t` the type of the basis,
  `n` the principle quantum number and `l` the angular momentum quantum number.
- `matrix.txt` for the input matrix elements on one particular k/R point.
  The sparse matrix-market format is adapted here.
  The indices of matrix is ordered as indicated by `basis_id.txt`, with azimuth quantum number from `-l` to `l`.

To run the program

```shell
abf_trans.exe lattice.txt basis_id.txt matrix.txt mode x1 x2 x3
```

mode should be either `R`/`K`, and `xN` are the components along direct/reciprocal lattice vectors.

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
