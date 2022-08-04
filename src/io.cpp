#include "io.h"
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>
// #include <iostream>

double decode_fraction(const string& frac_str)
{
    double frac = 0.0;
    const auto ite_slash = std::find(frac_str.cbegin(), frac_str.cend(), '/');
    if (ite_slash != frac_str.cend())
    {
        string s_nume, s_denu;
        s_nume = string(frac_str.begin(), ite_slash);
        s_denu = string(ite_slash+1, frac_str.end());
        // std::cout << s_nume << " " << s_denu << std::endl;
        frac = std::stod(s_nume) / std::stod(s_denu);
    }
    else
        frac = std::stod(frac_str);
    return frac;
}

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

matrix<cplxdb> read_mtx_cplxdb(const string &mtxfile)
{
    std::ifstream fin;
    fin.open(mtxfile);
    int nr, nc, nnz;
    matrix<cplxdb> cmat;
    if (fin)
    {
        string s1, s2, s3, s4;
        // skip two headlines
        getline(fin, s1);
        getline(fin, s1);
        fin >> s1 >> s2 >> s3;
        nr = std::stoi(s1);
        nc = std::stoi(s2);
        nnz = std::stoi(s3);
        cmat.resize(nr, nc);
        for (int i = 0; i < nnz; i++)
        {
            if (fin.eof())
                throw std::logic_error("EOF is reached unexpectedly in reading complex MTX");
            fin >> s1 >> s2 >> s3 >> s4;
            cmat(std::stoi(s1) - 1, std::stoi(s2) - 1) = cplxdb(std::stod(s3), std::stod(s4));
        }
    }
    return cmat;
}

matrix<double> read_mtx_double(const string &mtxfile)
{
    std::ifstream fin;
    fin.open(mtxfile);
    int nr, nc, nnz;
    matrix<double> mat;
    if (fin)
    {
        string s1, s2, s3;
        // skip two headlines
        getline(fin, s1);
        getline(fin, s1);
        fin >> s1 >> s2 >> s3;
        nr = std::stoi(s1);
        nc = std::stoi(s2);
        nnz = std::stoi(s3);
        mat.resize(nr, nc);
        for (int i = 0; i < nnz; i++)
        {
            if (fin.eof())
                throw std::logic_error("EOF is reached unexpectedly in reading complex MTX");
            fin >> s1 >> s2 >> s3;
            mat(std::stoi(s1) - 1, std::stoi(s2) - 1) = std::stod(s3);
        }
    }
    return mat;
}

void write_mtx_cplxdb(const matrix<cplxdb> &mat, const string &mtxfile,
                      const string &comment,
                      double threshold, bool row_first)
{
    std::ofstream fs;
    fs.open(mtxfile);
    int nr = mat.nr;
    int nc = mat.nc;
    int prec = 15;
    size_t nnz = 0;
    auto format = std::scientific;
    fs << "%%MatrixMarket matrix coordinate complex general" << std::endl;
    if (comment != "")
        fs << "% " << comment << std::endl;
    else
        fs << "%" << std::endl;
    // count non-zero values first
    for (int i = 0; i < mat.size(); i++)
    {
        auto v = mat.c[i];
        if ( fabs(v.real()) > threshold || fabs(v.imag()) > threshold )
            nnz++;
    }

    fs << nr << " " << nc << " " << nnz << std::endl;

    if (row_first)
    {
        for (int j = 0; j < nc; j++)
        {
            for (int i = 0; i < nr; i++)
            {
                auto v = mat(i, j);
                if ( fabs(v.real()) > threshold || fabs(v.imag()) > threshold )
                    fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v.real() << " " << std::showpoint << format << std::setprecision(prec) << v.imag() << "\n";
            }
        }
    }
    else
    {
        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nc; j++)
            {
                auto v = mat(i, j);
                if ( fabs(v.real()) > threshold || fabs(v.imag()) > threshold )
                    fs << i + 1 << " " << j + 1 << " " << std::showpoint << format << std::setprecision(prec) << v.real() << " " << std::showpoint << format << std::setprecision(prec) << v.imag() << "\n";
            }
        }
    }

}
