#include <utility>

//
// Created by David Winderl on 5/1/21.
//

#ifndef MITM_ALGORITHM_GATE_H
#define MITM_ALGORITHM_GATE_H


#include "utils.h"


using namespace std;



class Gate {

public:
    int register_position;
    int control_position;
    GateType type;

    Gate(int register_position, GateType type) {
        assert(type != CX);
        this->type = type;
        this->register_position = register_position;
        this->control_position = register_position;

    }

    Gate(int register_position, int control_position, GateType type) {
        assert(type == CX);
        this->type = type;
        this->register_position = register_position;
        this->control_position = control_position;
    }

    void get_unitary(xarrayc &out);

};


#endif //MITM_ALGORITHM_GATE_H
