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
    auto nabfs = get_number_of_abfs();
    n_tot_abfs = 0;
    for (const auto &nabf: nabfs)
    {
        start_index_atom.push_back(n_tot_abfs);
        n_tot_abfs += nabf;
    }
}

void ABF::get_abf_arlm(int abf_index, int &iat, int &irf, int &l, int &m) const
{
    if (abf_index >= n_tot_abfs)
        throw std::invalid_argument("requested index >= nbasis");
    for (int i = 1; i < atom_types.size(); i++)
    {
        if (abf_index > start_index_atom[i]) continue;
        iat = i - 1;
        int residual = abf_index - start_index_atom[iat];
        for (const auto &abf: map_type_abfs.at(atom_types[iat]))
        {
            residual -= abf.size();
            if (residual < 0)
            {
                irf = abf.irf;
                l = abf.l;
                m = - residual - abf.size() - l;
                break;
            }
        }
        break;
    }
}
