#include "gate.h"
#include "circuit.h"
#include "dispatch_queue.h"
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
    for (_GateType *tmp : set) {
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
    vector<_GateType> oneQubitGateSet = {h, s, sdg, t, tdg};
    vector<_GateType *> queue;
    _GateType tmp_archetype[size];
    generateArchetype(size, tmp_archetype);
    for (_GateType gate : oneQubitGateSet) {
        for (int j = 0; j < size; j++) {
            auto *new_vec = new _GateType[size];
            deepcopy(tmp_archetype, new_vec, size);
            new_vec[j] = gate;
            queue.push_back(new_vec);
        }
    }
    //queue.push_back(tmp_archetype);
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
        for (_GateType gate : oneQubitGateSet) {
            if (isGateInCirc(tmp, gate, size))
                continue;
            for (int j = 0; j < size; j++) {
                if (tmp[j] == id) {
                    auto *new_vec = new _GateType[size];
                    deepcopy(tmp, new_vec, size);
                    new_vec[j] = gate;
                    queue.push_back(new_vec);
                }
            }
        }
    }
    vector<Circuit> circ_out;
    circ_out.reserve(out.size());
    for (const auto &tmp : out) {
        auto tmp_circ = Circuit::fromDepthOneArray(tmp, size);
        tmp_circ.generateHash(static_vecs, m);
        circ_out.push_back(tmp_circ);
    }
    return circ_out;
}

void
generate_hashes(const vector<xarrayc> &static_vecs, int m, const xarrayc &u, const xarrayc &unitary_circ,
                xarrayc &hash,
                xarrayc &hash_inv) {
    xarrayc unitary = xt::linalg::dot(xt::transpose(xt::conj(unitary_circ)), u);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            auto tmp_inv = xt::linalg::dot(xt::linalg::dot(xt::transpose(xt::conj(static_vecs[i])), unitary),
                                           static_vecs[j]);
            auto tmp = xt::linalg::dot(xt::linalg::dot(xt::transpose(xt::conj(static_vecs[i])), unitary),
                                       static_vecs[j]);
            hash_inv(i, j) = xt::real(tmp_inv(0, 0));
            hash(i, j) = xt::real(tmp(0, 0));
        }
    }


}

void generateNextSet(const vector<Circuit> &allCombinations,
                     const vector<Circuit> &s_i,
                     vector<Circuit> &s_i_new,
                     vector<Circuit> &s_i_new_inv,
                     vector<xarrayc> &static_vecs,
                     int m, const xarrayc &unitary) {
    s_i_new.reserve(allCombinations.size() * s_i.size());
    s_i_new_inv.reserve(allCombinations.size() * s_i.size());
    for (const Circuit &v : allCombinations) {
        for (const Circuit &s : s_i) {
            Circuit circ = s.compose(v);
            Circuit circ1 = s.compose(v);
            xarrayc unitary =
            circ.generateHash(static_vecs, m);
            circ1.generateHashInv(static_vecs, m, unitary);
            s_i_new.emplace_back(circ);
            s_i_new_inv.emplace_back(circ1);
        }
    }
}

std::mutex mtx;

void generateNextSetAsync(const vector<Circuit> &allCombinations,
                          const vector<Circuit> &s_i,
                          vector<Circuit> &s_i_new,
                          vector<Circuit> &s_i_new_inv,
                          vector<xarrayc> &static_vecs,
                          int m, const xarrayc &unitary) {
    s_i_new.reserve(allCombinations.size() * s_i.size());
    s_i_new_inv.reserve(allCombinations.size() * s_i.size());

    dispatch_queue_t task_queue;
    for (const Circuit &v : allCombinations) {
        for (const Circuit &s : s_i) {
            task_queue.enqueue([&]() {
                Circuit circ = s.compose(v);
                Circuit circ1 = s.compose(v);
                circ.generateHash(static_vecs, m);
                circ1.generateHashInv(static_vecs, m, unitary);
                mtx.lock();
                s_i_new.emplace_back(circ);
                s_i_new_inv.emplace_back(circ1);
                mtx.unlock();
            });
        }
    }
    task_queue.wait();
}

bool hasIntersection(const vector<Circuit> &s, const vector<Circuit> &s_inv, const xarrayc &unitary,
                     Circuit &v,
                     Circuit &w) {
    vector<Circuit> intersection_s_inv, intersection_s;
    unsigned long min_size = INT_MAX;
    bool found_intersection = false;
    int i = 0, j = 0, s_size = s.size(), s1_size = s_inv.size();
    while (i < s_size && j < s1_size) {
        if (s[i] < s_inv[j]) {
            i++;
        } else if (s[i] > s_inv[j]) {
            j++;
        } else {
            auto w_ = s[i];
            auto v_ = s_inv[j];
            i++;
            j++;
            Circuit out = w_.compose(v_);
            if (out.gates.size() > min_size) {
                continue;
            } else {
                auto u_v = v_.getUnitary();
                u_v = xt::linalg::dot(xt::transpose(xt::conj(u_v)), unitary);
                auto u_w = w_.getUnitary();
                if (xt::allclose(u_v, u_w)) {
                    v = v_;
                    w = w_;
                    found_intersection = true;
                    min_size = out.gates.size();
                }
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


Circuit mitmAlgorithm(const xarrayc &unitary, int depth, int numQubits, int m = 2, bool verbose = true) {
    vector<xarrayc> static_vecs = fillStaticVecs(m, numQubits);
    vector<Circuit> s_i, s_i_inv_u;
    vector<Circuit> s_i_prev, s_i_inv_u_prev;
    const vector<Circuit> allCombinations = generateAllCombinations(numQubits, static_vecs, m);
    auto circ = Circuit(numQubits);
    s_i_prev.emplace_back(circ);
    s_i_inv_u_prev.emplace_back(circ);
    for (int i = 1; i < (int) depth / 2.0; ++i) {
        if (verbose)
            cout << "depth: " << i << endl;

        if (verbose)
            cout << "Generating all combinations" << endl;
        s_i.clear();
        s_i_inv_u.clear();
        if (depth < 2)
            generateNextSet(allCombinations, s_i_prev, s_i, s_i_inv_u, static_vecs, m, unitary);
        else
            generateNextSetAsync(allCombinations, s_i_prev, s_i, s_i_inv_u, static_vecs, m, unitary);
        if (verbose)
            cout << "Sorting..." << endl;
        std::sort(s_i.begin(), s_i.end());
        std::sort(s_i_inv_u.begin(), s_i_inv_u.end());
        if (verbose)
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
    if (verbose)
        cout << "Ran out of depth!" << endl;
    return Circuit(0);
}

int main() {
    xt::xarray<complex<double>> unitary_cy = {
            {1.0, 0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, -1.i},
            {0.0, 0.0, 1.i, 1.0}};
    clock_t begin = clock();
    auto circ = mitmAlgorithm(unitary_cy, 8, 2, 1, true);
    clock_t end = clock();
    cout << double(end - begin) / CLOCKS_PER_SEC << endl;
    cout << circ.toString() << endl;
    return 0;
}
