/*!
 @file spglib_wrapper.h
 @brief The file contains some wrapper functions of spglib utilities
 */
#pragma once
#include <spglib.h>
#include "matrix.h"
#include <vector>
#include <string>

using std::string;
using std::vector;

// a wrapper class for the original SpglibDataset struct
struct SpgDS_c
{
public:
    matrix<double> lattice;
    matrix<double> positions;
    vector<int> types;
    int spacegroup_number;
    int hall_number;
    string international_symbol;
    string hall_symbol;
    string pointgroup_symbol;
    string choice;
    matrix<double> transformation_matrix;
    vector<double> origin_shift;
    int n_operations;
    vector<matrix<int>> rotations;
    vector<vector<double>> translations;
    int n_atoms;
    vector<int> wyckoffs;
    vector<string> site_symmetry_symbols;
    vector<int> equivalent_atoms;
    vector<int> crystallographic_orbits;
    matrix<double> primitive_lattice;
    vector<int> mapping_to_primitive;
    int n_std_atoms;
    matrix<double> std_lattice;
    vector<int> std_types;
    matrix<double> std_positions;
    matrix<double> std_rotation_matrix;
    vector<int> std_mapping_to_primitive;

    SpgDS_c(const matrix<double> &latt_in, const matrix<double> &posi_frac_in,
            const vector<int> &types_in, const double symprec = 1.0e-5);
    ~SpgDS_c() {};
    // std_lattice is already the ideal one
    matrix<double> get_ideal_lattice() { return std_lattice; }
    const matrix<double> &get_ideal_lattice() const { return std_lattice; }

    void show() const;
};

SpglibDataset* wrapper_spg_get_dataset(const matrix<double> &latt,
                                       const matrix<double> &posi_frac,
                                       const vector<int> &types,
                                       const double symprec = 1.0e-5);

void show_spg_dataset(const SpglibDataset *dataset);

matrix<double> posi_frac_upon_symop(const SpglibDataset *dataset, const int &isym);
