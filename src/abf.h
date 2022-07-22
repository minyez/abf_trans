#pragma once
#include <map>
#include <vector>
#include "base.h"

struct abf_id
{
    int kind;
    int n;
    int l;
    abf_id(int kind_in, int n_in, int l_in): kind(kind_in), n(n_in), l(l_in) {};
    int size() const { return get_msize(l); }
};

std::vector<int> get_number_of_abfs(const std::map<int, std::vector<abf_id>> &map_type_abfs,
                                    const std::vector<int> &atypes);

int get_number_of_total_abfs(const std::map<int, std::vector<abf_id>> &map_type_abfs,
                             const std::vector<int> &atypes);

class ABF
{
private:
    std::vector<int> start_index_atom;
    int n_tot_abfs;
public:
    std::map<int, std::vector<abf_id>> map_type_abfs;
    std::vector<int> atom_types;
    ABF(const std::vector<int> &atom_types_in, const std::map<int, std::vector<abf_id>> &map_type_abfs_in);
    ~ABF() {};
    std::vector<int> get_number_of_abfs() const { return ::get_number_of_abfs(map_type_abfs, atom_types); }
    int get_number_of_total_abfs() const { return n_tot_abfs; }
    void get_abf_knlm(int abf_index, int &iat, int &kind, int &n, int &l, int &m) const;
    /* int get_abf_index(int iat, int kind, int n, int l, int m) const; */
};

extern std::map<int, std::vector<abf_id>> map_type_abfs;

