#include "io.h"
#include <fstream>
#include <stdexcept>

void read_cell(const string &cellfile, matrix<double> &out_latt,
               matrix<double> &out_posi_frac, vector<int> &out_types)
{
    if (out_latt.size() > 0 || out_posi_frac.size() > 0 || out_types.size() > 0)
        throw std::logic_error("cell has been read, please first clear it");

    std::ifstream fin;
    fin.open(cellfile);
    if (fin)
    {
        string s1, s2, s3, s4;
        out_latt.resize(3, 3);
        for (int i = 0; i < 3; i++)
        {
            fin >> s1 >> s2 >> s3;
            if (!fin)
                throw std::logic_error("error in reading the lattice");
            out_latt(i, 0) = std::stod(s1);
            out_latt(i, 1) = std::stod(s2);
            out_latt(i, 2) = std::stod(s3);
        }
        fin >> s1;
        int natoms = std::stoi(s1);
        out_posi_frac.resize(natoms, 3);
        for (int i = 0; i < natoms; i++)
        {
            fin >> s1 >> s2 >> s3 >> s4;
            if (!fin)
                throw std::logic_error("error in reading the atom positions and types");
            out_posi_frac(i, 0) = std::stod(s1);
            out_posi_frac(i, 1) = std::stod(s2);
            out_posi_frac(i, 2) = std::stod(s3);
            out_types.push_back(std::stoi(s4));
        }
    }
    else
        throw std::logic_error("fail to open cell input file, path: " + cellfile);
    fin.close();
}

void clean_cell(matrix<double> &in_latt, matrix<double> &in_posi_frac, vector<int> &in_types)
{
    in_latt.resize(0, 0);
    in_posi_frac.resize(0, 0);
    in_types.clear();
}

void read_abf_ids(const string &abffile, const set<int> &inequiv_types_verify, map<int, vector<abf_id>> &map_t_abf)
{
    if (!map_t_abf.empty())
        throw std::logic_error("atomic basis function map has been parsed");

    std::ifstream fin;
    fin.open(abffile);
    if(fin)
    {
        string t, l;
        while(!fin.eof())
        {
            fin >> t >> l;
            abf_id aid(std::stoi(l));
            map_t_abf[std::stoi(t)].push_back(aid);
        }
    }
    else
        throw std::logic_error("fail to open atomic basis input file, path: " + abffile);
    for (const auto &it: inequiv_types_verify)
    {
        if (map_t_abf.count(it) == 0)
            throw std::logic_error("missing atomic basis function information of atom type " + std::to_string(it));
    }
}
