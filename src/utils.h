#ifndef MITM_ALGORITHM_UTILS_H
#define MITM_ALGORITHM_UTILS_H


#include <vector>
#include <queue>
#include <math.h>
#include <iostream>
#include <future>


#include <xtensor/xarray.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor-blas/xlinalg.hpp>

enum GateType {
    H = 0,
    CX = 1,
    X = 2,
    Y = 3,
    Z = 4,
    S = 5,
    Sdg = 6,
    T = 7,
    Tdg = 8,
    I = 9 // note identity gate is not supported

};


/**
 * Internal Gatetype for parsing the control position of CNOT
 * at a depth 1 vector
 */
enum _GateType {
    h = 0,
    c = 1,
    x = 2,
    y = 3,
    z = 4,
    s = 5,
    sdg = 6,
    t = 7,
    tdg = 8,
    id = 9,
    _x = 10,

};


typedef xt::xarray<std::complex<double >> xarrayc;
typedef xt::xarray<double> xarrayf;




#endif //MITM_ALGORITHM_UTILS_H
