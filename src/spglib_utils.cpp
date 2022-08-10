#include "mathtools.h"
#include "spglib_utils.h"

SpgDS_c::SpgDS_c(const matrix<double> &latt_in, const matrix<double> &posi_frac_in,
                 const std::vector<int> &types_in, const double symprec_in) :
    lattice(latt_in), positions(posi_frac_in), types(types_in), symprec(symprec_in)
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
    n_std_atoms = dataset->n_std_atoms;

    std_rotation_matrix.resize(3, 3);
    transformation_matrix.resize(3, 3);
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            transformation_matrix(j, i) = dataset->transformation_matrix[i][j];
            std_rotation_matrix(j, i) = dataset->std_rotation_matrix[i][j];
        }
        origin_shift.push_back(dataset->origin_shift[i]);
    }

    rotations.resize(n_operations);
    translations.resize(n_operations);
    for (int iop = 0; iop < n_operations; iop++)
    {
        rotations[iop].resize(3, 3);
        translations[iop].resize(3);
        for (int i = 0; i < 3; i++)
        {
            translations[iop][i] = dataset->translations[iop][i];
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

    // tabulate inverse operations
    inverse_operation.resize(n_operations);
    for (int i = 0; i < n_operations; i++)
        inverse_operation[i] = -1;
    matrix<int> iden(3, 3);
    vec<double> zero(3);
    iden.set_diag(1);
    for (int i_op = 0; i_op < n_operations; i_op++)
    {
        if (inverse_operation[i_op] != -1) continue;
        for (int j_op = i_op; j_op < n_operations; j_op++)
        {
            if (rotations[i_op] * rotations[j_op] == iden)
            {
                vec<double> ff(3);
                for (int i = 0; i < 3; i++)
                {
                    // V_inv f + f_inv should be direct lattice vectors.
                    ff[i] = rotations[j_op](i, 0) * translations[i_op][0] +
                            rotations[j_op](i, 1) * translations[i_op][1] +
                            rotations[j_op](i, 2) * translations[i_op][2] + translations[j_op][i];
                    shift_to_unit(ff[i], 0.0, true);
                }
                if (ff == zero)
                {
                    inverse_operation[i_op] = j_op;
                    inverse_operation[j_op] = i_op;
                    break;
                }
            }
        }
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

void SpgDS_c::show(bool show_operations) const
{
    printf("International: %s (%d)\n", international_symbol.c_str(), spacegroup_number);
    printf("  Hall symbol: %s\n", hall_symbol.c_str());
    printf("  Point group: %s\n", pointgroup_symbol.c_str());
    printf("Transformation matrix:\n");
    for (int i = 0; i < 3; i++) {
      printf("  %12.6f %12.6f %12.6f\n", transformation_matrix(i, 0), transformation_matrix(i, 1), transformation_matrix(i, 2));
    }
    printf("Origin shift: [%f %f %f]^T\n", origin_shift[0],origin_shift[1],origin_shift[2]);
    printf("Symmetry operations (%d):\n", n_operations);
    if (show_operations)
        for (int i = 0; i < n_operations; i++)
            printf("%3d: %s (invop %3d) (%d)\n", i+1, get_operation_str(i).c_str(), inverse_operation[i]+1, is_proper(i));
    printf("Equivalent atoms:\n");
    for (int i = 0; i < n_atoms; i++) {
      printf("  %d -> %d (type %d)\n", i, equivalent_atoms[i], types[i]);
    }
}

void SpgDS_c::show_cell() const
{
    printf("Lattice parameter:\n");
    for (int i = 0; i < 3; i++) {
        printf("%f %f %f\n", lattice(i, 0), lattice(i, 1), lattice(i, 2));
    }
    printf("Atomic positions:\n");
    for (int i = 0; i < n_atoms; i++) {
        printf("%d: %f %f %f\n", types[i], positions(i, 0), positions(i, 1), positions(i, 2));
    }
}

string SpgDS_c::get_operation_str(int isymop) const
{
    char s[100];
    if (!(isymop < n_operations))
        throw std::invalid_argument("Required operation not exist");
    sprintf(s, "%2d %2d %2d %2d %2d %2d %2d %2d %2d [ %5.3f %5.3f %5.3f ]",
            rotations[isymop](0, 0), rotations[isymop](0, 1), rotations[isymop](0, 2),
            rotations[isymop](1, 0), rotations[isymop](1, 1), rotations[isymop](1, 2),
            rotations[isymop](2, 0), rotations[isymop](2, 1), rotations[isymop](2, 2),
            translations[isymop][0], translations[isymop][1], translations[isymop][2]
            );

    return string(s);
}

string SpgDS_c::get_operation_str_matform(int isymop) const
{
    char s[100];
    if (!(isymop < n_operations))
        throw std::invalid_argument("Required operation not exist");
    sprintf(s, "%2d %2d %2d | %5.3f\n%2d %2d %2d | %5.3f\n%2d %2d %2d | %5.3f",
            rotations[isymop](0, 0), rotations[isymop](0, 1), rotations[isymop](0, 2), translations[isymop][0],
            rotations[isymop](1, 0), rotations[isymop](1, 1), rotations[isymop](1, 2), translations[isymop][1],
            rotations[isymop](2, 0), rotations[isymop](2, 1), rotations[isymop](2, 2), translations[isymop][2]
            );

    return string(s);
}

bool SpgDS_c::is_proper(int isymop) const
{
    if (!(isymop < n_operations))
        throw std::invalid_argument("Required operation not exist");
    const matrix<int> &rot = rotations[isymop];
    int det = rot(0, 0) * (rot(1, 1) * rot(2, 2) - rot(1, 2) * rot(2, 1))
            + rot(0, 1) * (rot(1, 2) * rot(2, 0) - rot(1, 0) * rot(2, 2))
            + rot(0, 2) * (rot(1, 0) * rot(2, 1) - rot(1, 1) * rot(2, 0));
    if (det == 1)
        return true;
    if (det == -1)
        return false;
    throw std::logic_error("Determinant of rotation matrix is neither 1 or -1");
}
