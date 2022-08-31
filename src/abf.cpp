#include <set>
#include "abf.h"

std::map<int, std::vector<abf_id>> map_type_abfs;

std::vector<int> get_number_of_abfs(const std::map<int, std::vector<abf_id>> &map_type_abfs,
                                    const std::vector<int> &atypes)
{
    std::vector<int> nabfs;
    std::map<int, int> map_type_nabf;
    for (const auto &at: atypes)
    {
        int nabf = 0;
        if (map_type_nabf.count(at))
            nabf = map_type_nabf.at(at);
        else
        {
            if (0 == map_type_abfs.count(at))
                throw std::logic_error("Missing atom type in ABF mapping");
            for (const auto &abf: map_type_abfs.at(at))
            {
                nabf += abf.size();
            }
            map_type_nabf[at] = nabf;
        }
        nabfs.push_back(nabf);
    }
    return nabfs;
}

int get_number_of_total_abfs(const std::map<int, std::vector<abf_id>> &map_type_abfs,
                             const std::vector<int> &atypes)
{
    auto abfs = get_number_of_abfs(map_type_abfs, atypes);
    int n_tot_abfs = 0;
    for (const auto &i: abfs)
        n_tot_abfs += i;
    return n_tot_abfs;
}

ABF::ABF(const std::vector<int> &atom_types_in, const std::map<int, std::vector<abf_id>> &map_type_abfs_in)
        : atom_types(atom_types_in), map_type_abfs(map_type_abfs_in)
{
    /* // check duplicates */
    /* for (const auto &at_abfs: map_type_abfs_in) */
    /* { */
    /*     const auto & abfs = at_abfs.second; */
    /*     std::set<abf_id> set_abfs(abfs.cbegin(), abfs.cend()); */
    /*     if (abfs.size() != set_abfs.size()) */
    /*         throw std::invalid_argument("found duplicates in map_type_abfs input, please check"); */
    /* } */
    auto nabfs = get_number_of_abfs();
    n_tot_abfs = 0;
    for (const auto &nabf: nabfs)
    {
        start_index_atom.push_back(n_tot_abfs);
        n_tot_abfs += nabf;
        end_index_atom.push_back(n_tot_abfs - 1);
    }
}

void ABF::get_abf_arlm(int abf_index, int &iat, int &irf, int &l, int &m, const CODE_CHOICE &code) const
{
    if (abf_index >= n_tot_abfs)
        throw std::invalid_argument("requested index >= nbasis");
    for (int i = 0; i < atom_types.size(); i++)
        if (abf_index >= start_index_atom[i] && abf_index <= end_index_atom[i])
            iat = i;
    int residual = abf_index - start_index_atom[iat];
    for (int i = 0; i < map_type_abfs.at(atom_types[iat]).size(); i++)
    {
        auto abf = map_type_abfs.at(atom_types[iat])[i];
        residual -= abf.size();
        if (residual < 0)
        {
            residual += abf.size();
            irf = i;
            l = abf.l;
            switch (code)
            {
                case (CODE_CHOICE::ABACUS):
                    m = ((residual+1) / 2) * (residual%2? 1 : -1);
                    break;
                default:
                    m = residual - l;
            }
            break;
        }
    }
}

int ABF::get_abf_index(int iat, int irf, int m, const CODE_CHOICE &code) const
{
    if (!(iat < atom_types.size()))
        throw std::invalid_argument("atom index out of bound");
    int iabf = start_index_atom[iat];
    auto abfs = map_type_abfs.at(atom_types[iat]);
    if (!(irf < abfs.size()))
        throw std::invalid_argument("radfunc index out of bound");
    int l = abfs[irf].l;
    for (int i = 0; i < irf; i++)
        iabf += abfs[i].size();
    switch (code)
    {
        case (CODE_CHOICE::ABACUS):
            iabf += 2*std::abs(m) - int(m > 0);
            break;
        default:
            iabf += l + m;
    }
    return iabf;
}

void get_bloch_phase_convention(const CODE_CHOICE &cc, int &a, int &b)
{
    switch (cc)
    {
        case (CODE_CHOICE::ORIG):
            a = 1; b = 0; break;
        case (CODE_CHOICE::AIMS):
            a = -1; b = 0; break;
        case (CODE_CHOICE::ABACUS):
            a = 1; b = 0; break;
    }
}
