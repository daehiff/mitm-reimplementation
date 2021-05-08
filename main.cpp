#include "Gate.h"
#include "Circuit.h"

#include "utils.h"


using namespace std;
using namespace xt::linalg;

/**
 * Generate the archetype of an empty circuit
 * [i, i, i, ...].T
 * @param n
 * @param out
 */
void generateArchetype(int n, _GateType *out) {
    for (int j = 0; j < n; ++j) {
        out[j] = id;
    }
}

/**
 * Deepcopies an array of GateTypes_
 * @param gatevecin
 * @param gatevecout
 * @param size
 */
void deepcopy(const _GateType *gatevecin, _GateType *gatevecout, int size) {
    for (int i = 0; i < size; ++i) {
        gatevecout[i] = gatevecin[i];
    }
}

/**
 * Check if a certain gate is in the circuit
 * @param gatevec
 * @param gate
 * @param size
 * @return
 */
bool isGateInCirc(const _GateType *gatevec, _GateType gate, int size) {
    for (int i = 0; i < size; ++i) {
        if (gatevec[i] == gate)
            return true;
    }
    return false;
}

/**
 * Check if two Gatevecs are equal
 * @param gatevec
 * @param gatevec1
 * @param size
 * @return
 */
bool isEqual(const _GateType *gatevec, const _GateType *gatevec1, int size) {
    for (int i = 0; i < size; ++i) {
        if (gatevec[i] != gatevec1[i])
            return false;
    }
    return true;
}

bool inSet(const vector<_GateType *> &set, _GateType *circ, int size) {
    for (_GateType *tmp: set) {
        if (isEqual(tmp, circ, size)) {
            return true;
        }
    }
    return false;
}

/**
 * TODO make mem optimal
 * TODO paralellisation?
 * @param size
 * @return
 */
vector<Circuit> generateAllCombinations(int size, const vector<xarrayc> &static_vecs, int m) {
    vector<_GateType> oneQubitGateSet = {h, x, y, z, s, sdg, t, tdg};
    vector<_GateType *> queue;
    _GateType tmp_archetype[size];
    generateArchetype(size, tmp_archetype);
    queue.push_back(tmp_archetype);
    for (int j = 0; j < size; ++j) {
        for (int k = 0; k < size; ++k) {
            if (j != k) {
                auto *tmp = new _GateType[size];
                generateArchetype(size, tmp);
                tmp[j] = c;
                tmp[k] = _x;
                queue.push_back(tmp);
            }
        }
    }
    vector<_GateType *> out;
    while (!queue.empty()) {
        auto tmp = queue.back();
        queue.pop_back();
        if (!inSet(out, tmp, size)) {
            out.push_back(tmp);
        }
        for (_GateType gate: oneQubitGateSet) {
            if (isGateInCirc(tmp, gate, size))
                continue;
            for (int j = 0; j < size; j++) {
                if (tmp[j] == id) {
                    auto *new_vec = new _GateType[size]; // TODO allocate on stack
                    deepcopy(tmp, new_vec, size);
                    new_vec[j] = gate;
                    queue.push_back(new_vec);
                }
            }
        }
    }
    vector<Circuit> circ_out;
    circ_out.reserve(out.size());
    for (const auto &tmp: out) {
        auto tmp_circ = Circuit::fromDepthOneArray(tmp, size);
        tmp_circ.generateHash(static_vecs, m);
        circ_out.push_back(tmp_circ);
    }
    return circ_out;


}

// TODO make this multiprocessing to get more effective
void generateNextSet(const vector<Circuit> &allCombinations,
                     const vector<Circuit> &s_i,
                     vector<Circuit> &s_i_new,
                     vector<Circuit> &s_i_new_inv,
                     vector<xarrayc> &static_vecs,
                     int m, const xarrayc &unitary) {
    s_i_new.reserve(allCombinations.size() * s_i.size());
    s_i_new_inv.reserve(allCombinations.size() * s_i.size());
    for (const Circuit &v: allCombinations) {
        for (const Circuit &s: s_i) {
            Circuit circ = s.compose(v);
            Circuit circ1 = s.compose(v);
            circ.generateHash(static_vecs, m);
            circ1.generateHashInv(static_vecs, m, unitary);
            s_i_new.emplace_back(circ);
            s_i_new_inv.emplace_back(circ1);

        }
    }


}

bool hasIntersection(const vector<Circuit> &s, const vector<Circuit> &s_inv, const xarrayc &unitary,
                     Circuit &v,
                     Circuit &w) {
    vector<Circuit> intersection_s_inv, intersection_s;

    int i = 0, j = 0, s_size = s.size(), s1_size = s_inv.size();
    while (i < s_size && j < s1_size) {
        if (s[i] < s_inv[j]) {
            i++;
        } else if (s[i] > s_inv[j]) {
            j++;
        } else {
            intersection_s.push_back(s[i]);
            intersection_s_inv.push_back(s_inv[j]);
            i++;
            j++;
        }
    }
    int min_size = INT_MAX;
    bool found_intersection = false;
    // Assume, that the intersection is by far less due few key collisions
    for (auto w_: intersection_s) {
        for (auto v_: intersection_s_inv) {
            auto u_v = v_.getUnitary();
            u_v = xt::linalg::dot(xt::transpose(xt::conj(u_v)), unitary);

            auto u_w = w_.getUnitary();
            Circuit out = w_.compose(v_);
            if (xt::allclose(u_v, u_w) && out.gates.size() < min_size) {
                v = v_;
                w = w_;
                min_size = out.gates.size();
                found_intersection = true;
            }
        }
    }
    return found_intersection;
}


double fRand(double fMin, double fMax) {
    double f = (double) rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


vector<xarrayc> fillStaticVecs(int m, int size) {
    vector<xarrayc> static_vecs;
    static_vecs.reserve(m);
    for (int i = 0; i < m; ++i) {
        xarrayc tmp = xt::zeros<complex<double>>({(int) pow(2, size)});
        for (int j = 0; j < pow(2, size); ++j) {
            tmp(j) = 1. * fRand(-RAND_MAX, RAND_MAX) + 1.i * fRand(-RAND_MAX, RAND_MAX);
        }
        static_vecs.push_back(xt::eval(tmp));
    }
    return static_vecs;
}


Circuit mitmAlgorithm(const xarrayc &unitary, int depth, int numQubits, int m = 1) {
    vector<xarrayc> static_vecs = fillStaticVecs(m, numQubits);
    vector<Circuit> s_i, s_i_inv_u;
    vector<Circuit> s_i_prev, s_i_inv_u_prev;
    const vector<Circuit> allCombinations = generateAllCombinations(numQubits, static_vecs, m);
    auto circ = Circuit(numQubits);
    s_i_prev.emplace_back(circ);
    s_i_inv_u_prev.emplace_back(circ);
    for (int i = 1; i < (int) depth / 2.0; ++i) {
        cout << "depth: " << i << endl;

        cout << "Generating all combinations" << endl;
        s_i.clear();
        s_i_inv_u.clear();
        generateNextSet(allCombinations, s_i_prev, s_i, s_i_inv_u, static_vecs, m, unitary);
        cout << "Sorting..." << endl;
        std::sort(s_i.begin(), s_i.end());
        std::sort(s_i_inv_u.begin(), s_i_inv_u.end());


        cout << "Checking intersections" << endl;
        Circuit v = Circuit(numQubits);
        Circuit w = Circuit(numQubits);
        bool prev_intersection = hasIntersection(s_i, s_i_inv_u_prev, unitary, v, w);
        if (prev_intersection) {
            return w.compose(v);
        }

        bool current_intersection = hasIntersection(s_i, s_i_inv_u, unitary, v, w);
        if (current_intersection) {
            return w.compose(v);
        }
        s_i_prev = s_i;
        s_i_inv_u_prev = s_i_inv_u;

    }
    cout << "Ran out of depth!" << endl;
    return Circuit(0);
}

// TODO little sanaty checks, extends this as the alg is working :)
void circuit_test() {
    int numQubits = 3;
    Circuit circ = Circuit(numQubits);
    circ.addGate(CX, 1, 0);
    circ.addGate(H, 0);
    circ.addGate(X, 0);
    circ.addGate(Y, 1);
    circ.addGate(Z, 0);
    circ.addGate(T, 2);
    circ.addGate(Tdg, 1);
    circ.addGate(S, 2);
    circ.addGate(Sdg, 1);
    cout << circ.toString() << endl;
    auto u = circ.getUnitary();
    //cout << u << endl;
    cout << "Is Identity: " << endl;
    cout <<
         ((xt::allclose(xt::eye<complex<double>>((int) pow(2, numQubits)),
                        xt::linalg::dot(xt::transpose(xt::conj(u)), u)))
          ? "True" : "False")
         << endl;
}


void test_hash() {
    int m = 1, numQubits = 2;
    /*vector<xarrayc> static_vecs;
    static_vecs.push_back({-23.995845 + 38.832744i, 66.930692 - 93.855779i,
                           -8.07673 + 5.400387i, 34.298768 - 98.603628i});

*/
    vector<xarrayc> static_vecs = fillStaticVecs(m, numQubits);
    xarrayc unitary = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, -1.i},
            {0.0, 0.0, 1.i, 0.0}};

    Circuit circ = Circuit(2);
    circ.addGate(CX, 1, 0);
    circ.addGate(S, 1);
    circ.generateHashInv(static_vecs, m, unitary);

    Circuit circ1 = Circuit(2);
    circ1.addGate(Sdg, 1);
    circ1.generateHash(static_vecs, m);

    auto c1_u = circ1.getUnitary();
    auto c_u = circ.getUnitary();
    c_u = xt::linalg::dot(xt::transpose(xt::conj(c_u)), unitary);
    cout << (xt::real(circ1.hash) < xt::real(circ.hash)) << endl;
    cout << (xt::real(circ1.hash) > xt::real(circ.hash)) << endl;
    cout << ((xt::allclose(c1_u, c_u)) ? "True" : "False") << endl;
}


void test_set_generation() {
    int m = 1;
    int numQubits = 2;
    vector<xarrayc> static_vecs = fillStaticVecs(m, numQubits);
    xarrayc unitary = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, -1.i},
            {0.0, 0.0, 1.i, 0.0}};

    vector<Circuit> s1;
    vector<Circuit> s2;
    Circuit circ = Circuit(2);
    circ.addGate(CX, 1, 0);
    circ.addGate(S, 1);
    circ.generateHashInv(static_vecs, m, unitary);
    s1.push_back(circ);
    Circuit circ1 = Circuit(2);
    circ1.addGate(Sdg, 1);
    circ1.generateHash(static_vecs, m);
    s2.push_back(circ1);
    Circuit v = Circuit(numQubits);
    Circuit w = Circuit(numQubits);
    auto has_intersection = hasIntersection(s2, s1, unitary, v, w);
    cout << ((has_intersection) ? "True" : "False") << endl;
}


void test_mitm_algorithm() {
    // unitary w.r.t CY => s*, CX, s
    xt::xarray<complex<double>> unitary_cy = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, -1.i},
            {0.0, 0.0, 1.i, 0.0}};
    auto circ = mitmAlgorithm(unitary_cy, 10, 2, 1);
    cout << circ.toString() << endl;
    assert(circ.gates.size() == 3);
    assert(circ.gates[0].type == Sdg && circ.gates[0].register_position == 1);
    assert(circ.gates[1].type == CX && circ.gates[1].register_position == 1);
    assert(circ.gates[2].type == S && circ.gates[2].register_position == 1);
    // TODO find other test-cases from paper and test them against circuit
}

int main() {
    test_mitm_algorithm();
    return 0;
}
