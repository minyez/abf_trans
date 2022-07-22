#pragma once
#include <map>
#include <vector>
#include "base.h"

// ABF object for generalization
struct abf_id
{
private:
    int uid; // unique number for identity comparison
public:
    int l;
    abf_id(int l_in): l(l_in) { uid = l; }
    int size() const { return get_msize(l); }
    bool operator==(const abf_id &a2) const { return uid == a2.uid; }
    bool operator<(const abf_id &a2) const { return uid < a2.uid; }
    bool operator<=(const abf_id &a2) const { return uid <= a2.uid; }
    bool operator>(const abf_id &a2) const { return uid > a2.uid; }
    bool operator>=(const abf_id &a2) const { return uid >= a2.uid; }
    bool operator!=(const abf_id &a2) const { return !(*this == a2); }
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
    void get_abf_arlm(int abf_index, int &iat, int &irf, int &l, int &m) const;
    /* int get_abf_index(int iat, int irf, int l, int m) const; */
};

extern std::map<int, std::vector<abf_id>> map_type_abfs;

