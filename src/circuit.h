//
// Created by David Winderl on 5/1/21.
//

#ifndef MITM_ALGORITHM_CIRCUIT_H
#define MITM_ALGORITHM_CIRCUIT_H

#include "utils.h"
#include "circuit.h"
#include "gate.h"

using namespace std;

class Circuit {
public:
    vector<string> registerNames;
    vector<Gate> gates;
    int numQubits;
    xarrayf hash;

    Circuit(int numQubits) {
        this->numQubits = numQubits;
        this->gates = vector<Gate>();
        this->registerNames = vector<string>();
    }

    bool operator<(const Circuit &other) const {
        return xt::all(this->hash < other.hash);
    }

    bool operator>(const Circuit &other) const {
        return xt::all(this->hash > other.hash);
    }


    void addGate(GateType type, int q_register, int c_register);

    vector<GateType> getGateSet();

    xarrayc getUnitary() const;

    void toVectorRepresentation(vector<_GateType> &out);

    void generateHash(const vector<xarrayc> &static_vecs, int m);


    void generateHashInv(const vector<xarrayc> &static_vecs, int m, const xarrayc &unitary);

    Circuit inverse();

    string toString();

    static Circuit fromDepthOneArray(_GateType *depthOneVector, int size);

    bool isEqual(Circuit other);

    Circuit compose(const Circuit &other) const;

    void addGate(GateType type, int q_register);

    string toqcFormat();
};


#endif //MITM_ALGORITHM_CIRCUIT_H
