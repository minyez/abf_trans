#include "spglib_utils.h"

SpgDS_c::SpgDS_c(const matrix<double> &latt_in, const matrix<double> &posi_frac_in,
                 const std::vector<int> &types_in, const double symprec) : 
    lattice(latt_in), positions(posi_frac_in), types(types_in)
{
    SpglibDataset *dataset = wrapper_spg_get_dataset(latt_in, posi_frac_in, types_in, symprec);

    spacegroup_number = dataset->spacegroup_number;
    hall_number = dataset->hall_number;
    international_symbol = dataset->international_symbol;
    hall_symbol = dataset->hall_symbol;
    pointgroup_symbol = dataset->pointgroup_symbol;
    choice = dataset->choice;
    n_operations = dataset->n_operations;
    n_atoms = dataset->n_atoms;

    std_rotation_matrix.resize(3, 3);
    transformation_matrix.resize(3, 3);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            transformation_matrix(i, j) = dataset->transformation_matrix[i][j];
            std_rotation_matrix(i, j) = dataset->std_rotation_matrix[i][j];
        }
        origin_shift.push_back(dataset->origin_shift[i]);
    }

    rotations.resize(n_operations);
    translations.resize(n_operations);
    for (int iop = 0; iop < n_operations; iop++)
    {
        rotations[iop].resize(3, 3);
        for (int i = 0; i < 3; i++)
        {
            translations[iop].push_back(dataset->translations[iop][i]);
            for (int j = 0; j < 3; j++)
                rotations[iop](i, j) = dataset->rotations[iop][i][j];
        }
    }

    for (int ia = 0; ia < n_atoms; ia++)
    {
        wyckoffs.push_back(dataset->wyckoffs[ia]);
        site_symmetry_symbols.push_back(dataset->site_symmetry_symbols[ia]);
        equivalent_atoms.push_back(dataset->equivalent_atoms[ia]);
        crystallographic_orbits.push_back(dataset->crystallographic_orbits[ia]);
        mapping_to_primitive.push_back(dataset->mapping_to_primitive[ia]);
    }

    primitive_lattice.resize(3, 3);
    std_lattice.resize(3, 3);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            primitive_lattice(i, j) = dataset->primitive_lattice[j][i];
            std_lattice(i, j) = dataset->std_lattice[j][i];
        }

    std_positions.resize(n_std_atoms, 3);
    for (int ia = 0; ia < n_std_atoms; ia++)
    {
        for (int i = 0; i < 3; i++)
            std_positions(ia, i) = dataset->std_positions[ia][i];
        std_types.push_back(dataset->std_types[ia]);
        std_mapping_to_primitive.push_back(dataset->std_mapping_to_primitive[ia]);
    }

    spg_free_dataset(dataset);
}

SpglibDataset* wrapper_spg_get_dataset(const matrix<double> &latt_in,
                                       const matrix<double> &posi_frac_in,
                                       const vector<int> &types_in,
                                       const double symprec)
{
    SpglibDataset *dataset;
    const int natoms = types_in.size();
    int arr_types[natoms];
    double arr_latt[3][3];
    // spglib require basis vector in column format
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            arr_latt[i][j] = latt_in(j, i);
    double arr_posi_frac[natoms][3];
    for (int ia = 0; ia < natoms; ia++)
    {
        for (int i = 0; i < 3; i++)
            arr_posi_frac[ia][i] = posi_frac_in(ia, i);
        arr_types[ia] = types_in[ia];
    }
    dataset = spg_get_dataset(arr_latt, arr_posi_frac, arr_types, natoms, symprec);
    return dataset;
}

void SpgDS_c::show() const
{
    printf("International: %s (%d)\n", international_symbol.c_str(), spacegroup_number);
    printf("  Hall symbol: %s\n", hall_symbol.c_str());
    printf("  Point group: %s\n", pointgroup_symbol.c_str());
    printf("\n");
    printf("Symmetry operations (%d):\n", n_operations);
    for (int i = 0; i < n_operations; i++)
    {
        printf("%3d: %2d %2d %2d %2d %2d %2d %2d %2d %2d [ %f %f %f ]\n", i+1,
               rotations[i](0, 0), rotations[i](0, 1), rotations[i](0, 2),
               rotations[i](1, 0), rotations[i](1, 1), rotations[i](1, 2),
               rotations[i](2, 0), rotations[i](2, 1), rotations[i](2, 2),
               translations[i][0], translations[i][1], translations[i][2]
               );
    }
    printf("Equivalent atoms:\n");
    for (int i = 0; i < n_atoms; i++) {
      printf("  %d ->%d (type %d)\n", i, equivalent_atoms[i], types[i]);
    }
}


// routine adapted from https://github.com/spglib/spglib/blob/develop/example/example.c#L880
void show_spg_dataset(const SpglibDataset *dataset)
{
    int i, j, size;
    char ptsymbol[6];
    int pt_trans_mat[3][3];
    const char *wl = "abcdefghijklmnopqrstuvwxyz";

    printf("International: %s (%d)\n", dataset->international_symbol, dataset->spacegroup_number);
    printf("  Hall symbol: %s\n", dataset->hall_symbol);
    spg_get_pointgroup(ptsymbol, pt_trans_mat, dataset->rotations,
                       dataset->n_operations);
    printf("  Point group: %s\n", ptsymbol);
    // printf("Transformation matrix:\n");
    // for (i = 0; i < 3; i++) {
    //   printf("%f %f %f\n", dataset->transformation_matrix[i][0],
    //          dataset->transformation_matrix[i][1],
    //          dataset->transformation_matrix[i][2]);
    // }
    // printf("Origin shift: [%f %f %f]^T\n",
    //        dataset->origin_shift[0],dataset->origin_shift[1],dataset->origin_shift[2]);
    // printf("Wyckoff letters:\n");
    // for (i = 0; i < dataset->n_atoms; i++) {
    //   printf("%c ", wl[dataset->wyckoffs[i]]);
    // }
    // printf("\n");
    // printf("Standard atoms: %d. Poistions\n", dataset->n_std_atoms);
    // for (i = 0; i < dataset->n_std_atoms; i++) {
    //   printf("%d: %f %f %f\n", i+1, dataset->std_positions[i][0],  dataset->std_positions[i][1], dataset->std_positions[i][2]);
    // }
    // printf("\n");
    // printf("Equivalent atoms:\n");
    // for (i = 0; i < dataset->n_atoms; i++) {
    //   printf("%d ", dataset->equivalent_atoms[i]);
    // }
    // printf("\n");

    for (i = 0; i < dataset->n_operations; i++)
    {
        printf("%3d: %2d %2d %2d %2d %2d %2d %2d %2d %2d [ %f %f %f ]\n", i+1,
               dataset->rotations[i][0][0], dataset->rotations[i][0][1], dataset->rotations[i][0][2],
               dataset->rotations[i][1][0], dataset->rotations[i][1][1], dataset->rotations[i][1][2],
               dataset->rotations[i][2][0], dataset->rotations[i][2][1], dataset->rotations[i][2][2],
               dataset->translations[i][0], dataset->translations[i][1], dataset->translations[i][2]
               );
      // printf("--- %d ---\n", i + 1);
      // for (j = 0; j < 3; j++) {
      //   printf("%2d %2d %2d\n", dataset->rotations[i][j][0],
      //          dataset->rotations[i][j][1], dataset->rotations[i][j][2]);
      // }
      // printf("%f %f %f\n", dataset->translations[i][0],
      //        dataset->translations[i][1], dataset->translations[i][2]);
    }
}
