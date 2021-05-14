#include "pybind11/pybind11.h"

#include "xtensor/xmath.hpp"
#include "xtensor/xarray.hpp"

#define FORCE_IMPORT_ARRAY

#include "xtensor-python/pyarray.hpp"
#include "xtensor-python/pyvectorize.hpp"

#include <iostream>
#include <numeric>
#include <cmath>
#include "mitm.cpp"

#include "gate.cpp"
#include "circuit.cpp"

namespace py = pybind11;

std::string mitm_algorithm(xt::pyarray<std::complex<double>> &unitary, int max_depth, int m = 1, bool verbose = true) {
    if (unitary.shape()[0] != unitary.shape()[1])
        throw std::runtime_error("Requires a unitary of equal shape");
    int numQubits = (int) std::log2(unitary.shape()[0]);
    if (std::pow(2, numQubits) != unitary.shape()[0])
        throw std::runtime_error("Provided unitary has no shape of 2**n");
    auto circ = mitmAlgorithm(unitary, max_depth, numQubits, m, verbose);
    if (circ.gates.size() == 0)
        return "";
    return circ.toqcFormat();
}

PYBIND11_MODULE(mitm, m) {
    xt::import_numpy();
    m.doc() = R"pbdoc(
        An example xtensor extension

        .. currentmodule:: mitm

        .. autosummary::
           :toctree: _generate

           
    )pbdoc";

    m.def("mitm_algorithm", mitm_algorithm, "Mitm Algorithm as a c++ implementation");
}
