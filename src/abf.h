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
    int size() { return get_msize(l); }
};

extern std::map<int, std::vector<abf_id>> map_type_abfs;

