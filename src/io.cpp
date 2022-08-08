#include "io.h"
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

CODE_CHOICE code_choice = CODE_CHOICE::ORIG;
KRMODE krmode = KRMODE::K;

CODE_CHOICE parse_code_choice(const string &cc_in)
{
    if (cc_in == "orig" || cc_in == "original")
        return CODE_CHOICE::ORIG;
    else if (cc_in == "aims" || cc_in == "fhiaims")
        return CODE_CHOICE::AIMS;
    throw std::invalid_argument("Unknown code choice: " + cc_in);
}

KRMODE parse_krmode(const string &mode_in)
{
    if (mode_in == "r" || mode_in == "R")
        return KRMODE::R;
    else if (mode_in == "k" || mode_in == "K")
        return KRMODE::K;
    throw std::invalid_argument("Unknown mode: " + mode_in);
}

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
        string t, l, n;
        while(!fin.eof())
        {
            fin >> t >> l >> n;
            if (fin.eof())
                break;
            // std::cout << t << l << n << std::endl;
            abf_id aid(std::stoi(l));
            for (int i = 0; i < std::stoi(n); i++)
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

void read_matrix_inputs(const string &mat_inputs_fn, std::array<int, 3> &ngs, vector<vec<double>> &krpoints, vector<string> &mtxfns, vector<matrix<cplxdb>> &matrices)
{
    if (!krpoints.empty() || !mtxfns.empty() || !matrices.empty())
        throw std::logic_error("matrix inputs already parsed, please clean up first");

    std::ifstream fin;
    fin.open(mat_inputs_fn);
    if(fin)
    {
        int n = 0;
        string s1, s2, s3, s4, s5;
        fin >> s1 >> s2 >> s3 >> s4 >> s5;
        code_choice = parse_code_choice(s1);
        krmode = parse_krmode(s2);
        ngs[0] = stoi(s3);
        ngs[1] = stoi(s4);
        ngs[2] = stoi(s5);
        while (!fin.eof())
        {
            fin >> s1 >> s2 >> s3 >> s4;
            if (fin.eof()) break;
            n++;
            double varr[3] { decode_fraction(s1), decode_fraction(s2), decode_fraction(s3) };
            vec<double> v(3, varr);
            krpoints.push_back(v);
            mtxfns.push_back(s4);
            if (krmode == KRMODE::K)
                printf("Reading %d-th matrix at k-point (%6.3f, %6.3f, %6.3f) from file: %s\n", n, v[0], v[1], v[2], s4.c_str());
            if (krmode == KRMODE::R)
                printf("Reading %d-th matrix at R-point (%4.1f, %4.1f, %4.1f) from file: %s\n", n, v[0], v[1], v[2], s4.c_str());
            matrices.push_back(read_mtx_cplxdb(s4));
        }
        printf("%d files read\n", n);
    }
}

void clear_matrix_inputs(vector<vec<double>> &krpoints, vector<string> &matfns,
                         vector<matrix<cplxdb>> &matrices)
{
    krpoints.clear();
    matfns.clear();
    matrices.clear();
}
