//
// Created by David Winderl on 5/1/21.
//

#include "gate.h"

#include <iostream>
#include <xtensor.hpp>


using namespace std;
using namespace xt;

void Gate::get_unitary(xarrayc &out) {
            complex<double> tdg, t;
    switch (this->type) {
        case H:
            out = 1.0 / sqrt(2) * xarray<double>({{1.0, 1.0},
                                                  {1.0, -1.0}});
            break;
        case CX:
            out = xarray<double>({{1.0, 0.0, 0.0, 0.0},
                                  {0.0, 1.0, 0.0, 0.0},
                                  {0.0, 0.0, 0.0, 1.0},
                                  {0.0, 0.0, 1.0, 0.0}});
            break;
        case X:
            out = xarray<double>({{0.0, 1.0},
                                  {1.0, 0.0}});
            break;

        case Y:
            out = xarray<complex<double>>({{0.0, 1.0i},
                                           {-1i, 0.0}});
            break;
        case Z:
            out = xarray<double>({{1.0, 0.0},
                                  {0.0, -1.0}});
            break;
        case S:
            out = xarray<complex<double>>({{1.0, 0.0},
                                           {0.0, 1.0i}});
            break;
        case Sdg:
            out = xarray<complex<double>>({{1.0, 0.0},
                                           {0.0, -1.0i}});
            break;
        case T:
            t = std::cos(M_PI / 4.0) + std::sin(M_PI / 4.0) * 1.0i;
            out = xarray<complex<double>>({{1.0, 0.0},
                                           {0.0, t}});
            break;
        case Tdg:
            tdg = std::cos(M_PI / 4.0) - std::sin(M_PI / 4.0) * 1.0i;
            out = xarray<complex<double>>({{1.0, 0.0},
                                           {0.0, tdg}});
            break;
        case I:
            out = xt::eye(2);
            break;

        default:
            throw std::logic_error("Unknown gate type");

    }
}

