#include "matrix.h"
#include <array>

// the kgrids are always Gamma-centered
matrix<double> get_kgrids(std::array<int, 3> nks);
