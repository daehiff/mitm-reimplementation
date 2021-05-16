#include "circuit.h"

//
// Created by David Winderl on 5/1/21.
//


using namespace std;


std::string gateToString(GateType gate) {
    switch (gate) {

        case H:
            return "h ";
        case CX:
            return "cx";
        case X:
            return "x ";
        case Y:
            return "y ";
        case Z:
            return "z ";
        case S:
            return "s ";
        case Sdg:
            return "s*";
        case T:
            return "t ";
        case Tdg:
            return "t*";
        case I:
            return "i ";
    }
}

/**
 * TODO kron only needs to get executed once => should save performance (validate) and method gets simple
 * @param gate
 * @param circuitSize
 * @return
 */
void computeNFoldUnitary(Gate gate, int circuitSize, xarrayc &out) {
    int n = (gate.type == CX) ? circuitSize - 2 : circuitSize - 1;
    gate.get_unitary(out);
    if (n == 0)
        return;
    xarrayc identity = xt::eye((int) std::pow(2, n));
    out = xt::linalg::kron(out, identity);
}

vector<int> getOrder(Gate gate, int n) {
    vector<int> out;
    for (int i = 0; i < n; ++i) {
        out.push_back(i);
    }
    if (gate.type == CX) {
        int tmp = out[gate.control_position];
        int tmp1 = out[gate.register_position];
        out[gate.control_position] = 0;
        out[0] = tmp;
        out[gate.register_position] = 1;
        out[1] = tmp1;
    } else {
        int tmp = out[gate.register_position];
        out[gate.register_position] = 0;
        out[0] = tmp;
    }
    return out;

}


void reorderGates(xarrayc &gate, const vector<int> &permutation) {
    int n = permutation.size();
    vector<int> tmp;
    tmp.reserve(permutation.size() * 2);
    for (auto i : permutation) {
        tmp.push_back(i);
    }
    for (auto i : permutation) {
        tmp.push_back(i + n);
    }
    vector<int> twos;
    twos.reserve(2 * n);
    for (int i = 0; i < 2 * n; ++i) {
        twos.push_back(2);
    }
    int pow_ = pow(2, n);
    gate = transpose(gate.reshape(twos), tmp);
    gate = gate.reshape({pow_, pow_});
}


xarrayc Circuit::getUnitary() const {
    xarrayc out = xt::eye(pow(2, this->numQubits));
    xarrayc gate_unitary;
    for (auto gate: this->gates) {
        computeNFoldUnitary(gate, this->numQubits, gate_unitary);
        vector<int> order = getOrder(gate, this->numQubits);
        reorderGates(gate_unitary, order);
        out = xt::linalg::dot(gate_unitary, out);
    }
    return out;
}

void Circuit::addGate(GateType type, int q_register, int c_register) {
    if (type != CX) {
        auto gate = Gate(q_register, type);
        this->gates.push_back(gate);
    } else {
        auto gate = Gate(q_register, c_register, type);
        this->gates.push_back(gate);
    }

}

void Circuit::addGate(GateType type, int q_register) {
    this->addGate(type, q_register, -1);
}

vector<GateType> Circuit::getGateSet() {
    auto gateSet = vector<GateType>();
    for (auto gate: this->gates) {
        if (std::find(gateSet.begin(), gateSet.end(), gate.type) == gateSet.end())
            gateSet.push_back(gate.type);
    }
    return gateSet;
}

int getIndex(_GateType *v, _GateType gate, int size) {
    for (int i = 0; i < size; ++i) {
        if (v[i] == gate)
            return i;
    }
    return -1;
}

GateType convertType(_GateType type) {
    switch (type) {

        case h:
            return H;
        case c:
            throw std::invalid_argument("Would assume this cannot be used here");
        case x:
            return X;
        case y:
            return Y;
        case z:
            return Z;
        case s:
            return S;
        case sdg:
            return Sdg;
        case t:
            return T;
        case tdg:
            return Tdg;
        case id:
            return I;
        case _x:
            throw std::invalid_argument("Would assume this cannot be used here");
    }
}

Circuit Circuit::fromDepthOneArray(_GateType *depthOneVector, int size) {
    Circuit out = Circuit(size);
    vector<_GateType> oneQubitGateSet = {h, x, y, z, s, sdg, t, tdg};
    for (auto gate: oneQubitGateSet) {
        int index = getIndex(depthOneVector, gate, size);
        if (index != -1) {
            GateType type = convertType(gate);
            out.addGate(type, index);
        }
    }
    int indexC = getIndex(depthOneVector, c, size);
    int indexX = getIndex(depthOneVector, _x, size);

    if (indexC != -1) {
        assert(indexX != -1);
        out.addGate(CX, indexX, indexC);
    }
    return out;
}


void Circuit::toVectorRepresentation(vector<_GateType> &out) {
    for (auto gate :this->gates) {
        if (gate.type == CX) {
            out[gate.control_position] = c;
            out[gate.register_position] = _x;
        } else
            out[gate.register_position] = static_cast<_GateType >(gate.type);
    }

}


string Circuit::toString() {
    vector<string> lines;
    for (int j = 0; j < this->numQubits; ++j) {
        lines.push_back("q_" + to_string(j) + ": ");
    }

    for (auto gate: this->gates) {
        for (int j = 0; j < this->numQubits; ++j) {
            if (gate.type == CX) {
                if (j == gate.control_position)
                    lines[j] += " *  ";
                else if (j == gate.register_position)
                    lines[j] += " " + gateToString(gate.type) + " ";
                else
                    lines[j] += " i  ";
            } else {
                if (j == gate.control_position)
                    lines[j] += " " + gateToString(gate.type) + " ";
                else
                    lines[j] += " i  ";
            }
        }
    }

    string out;
    for (auto line : lines) {
        out += line + "\n";
    }
    return out;

}

bool Circuit::isEqual(Circuit other) {
    xarrayc thisUni = this->getUnitary();
    xarrayc otheUni = other.getUnitary();

    return xt::allclose(thisUni, otheUni);
}

Circuit Circuit::compose(const Circuit &other) const {
    Circuit out(this->numQubits);
    out.gates.reserve(this->gates.size() + other.gates.size());
    for (Gate gate: this->gates)
        out.addGate(gate.type, gate.register_position, gate.control_position);
    for (Gate gate: other.gates)
        out.addGate(gate.type, gate.register_position, gate.control_position);
    return out;
}

void Circuit::generateHash(const vector<xarrayc> &static_vecs, int m) {
    xarrayc unitary = this->getUnitary();
    this->hash = xt::zeros<double>({m, m});
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            auto tmp = xt::linalg::dot(xt::linalg::dot(xt::transpose(xt::conj(static_vecs[i])), unitary),
                                       static_vecs[j]);
            this->hash(i, j) = xt::real(tmp(0, 0));
        }
    }
}


void Circuit::generateHashInv(const vector<xarrayc> &static_vecs, int m, const xarrayc &u) {
    xarrayc unitary_circ = this->getUnitary();
    xarrayc unitary = xt::linalg::dot(xt::transpose(xt::conj(unitary_circ)), u);
    this->hash = xt::zeros<double>({m, m});
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            auto tmp = xt::linalg::dot(xt::linalg::dot(xt::transpose(xt::conj(static_vecs[i])), unitary),
                                       static_vecs[j]);
            this->hash(i, j) = xt::real(tmp(0, 0));
        }
    }
}

Circuit Circuit::inverse() {
    Circuit out = Circuit(this->numQubits);
    for (auto it = this->gates.rbegin(); it != this->gates.rend(); ++it) {
        switch (it->type) {
            case H:
            case CX:
            case X:
            case Y:
            case Z:
            case I:
                out.addGate(it->type, it->register_position, it->control_position);
                break;
            case S:
                out.addGate(Sdg, it->register_position, it->control_position);
                break;
            case Sdg:
                out.addGate(S, it->register_position, it->control_position);
                break;
            case T:
                out.addGate(Tdg, it->register_position, it->control_position);
                break;
            case Tdg:
                out.addGate(T, it->register_position, it->control_position);
                break;
            default:
                throw std::logic_error("Unknown gate type");
        }
    }
    return out;
}


std::string gateToQCString(GateType gate) {
    switch (gate) {

        case H:
            return "H ";
        case CX:
            return "tof";
        case X:
            return "X";
        case Y:
            return "Y";
        case Z:
            return "Z";
        case S:
            return "S";
        case Sdg:
            return "S*";
        case T:
            return "T";
        case Tdg:
            return "T*";
        case I:
            return "";
    }
}

string Circuit::toqcFormat() {
    string out;
    out += ".v ";
    for (int i = 0; i < this->numQubits; i++) {
        out += ("q" + to_string(i) + " ");
    }
    out += "\n";

    out += ".i ";
    for (int i = 0; i < this->numQubits; i++) {
        out += ("q" + to_string(i) + " ");
    }
    out += "\n";
    out += "\n";
    out += "BEGIN\n";
    for (auto gate: this->gates) {
        string qubit_position;
        if (gate.register_position != gate.control_position)
            qubit_position = "q" + to_string(gate.register_position) + " q" + to_string(gate.control_position);
        else
            qubit_position = "q" + to_string(gate.register_position);

        string gate_name = gateToQCString(gate.type);
        out += (gate_name + " " + qubit_position + "\n");

    }
    out += "END\n";
    return out;
}





