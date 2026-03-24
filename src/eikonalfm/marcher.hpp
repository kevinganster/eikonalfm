#pragma once

#include <limits>
#include <stddef.h>

// memory leak detection
// #ifdef _DEBUG
//	#include "debugtests.h"
// #endif

// "sizes" used (signed and unsigned)
using ssize = long;
using usize = unsigned long;

constexpr double INF = std::numeric_limits<double>::infinity();
constexpr double N_INF = std::numeric_limits<double>::lowest();
constexpr char UNKNOWN = 0;
constexpr char KNOWN = 1;
constexpr char FRONT = 2;

class MarcherInfo {
public:
    // number of dimensions
    int const ndim;
    // shape (of the mesh)
    usize const* const shape;
    // total number of grid points
    usize size;

    MarcherInfo(int const ndim, usize const* shape);

    virtual ~MarcherInfo() {}

    // No-Op functions for sensitivities in the base class
    virtual void store_sequence(const usize x) {}
    virtual void store_order(const int dim, const usize x, const char o) {}
    virtual auto get_sequence() -> usize* { return nullptr; }
    virtual auto get_order() -> char* { return nullptr; }
};

class SensitivityInfo: public MarcherInfo {
public:
    // execution sequence
    usize counter = 0;
    usize* sequence;
    // finite difference order
    char* orders;

    SensitivityInfo(int const ndim, usize const* shape): MarcherInfo{ndim, shape} {
        sequence = new usize[size];
        orders = new char[ndim * size];
    }

    ~SensitivityInfo() override {
        // the pointers are deleted in cfm.cpp
        // delete[] sequence;
        // delete[] orders;
    }

    void store_sequence(usize const x) override;
    void store_order(int const dim, usize const x, char const o) override;

    auto get_sequence() -> usize* override { return sequence; }
    auto get_order() -> char* override { return orders; }
};


class Marcher {
protected:
    // slowness field
    double const* const c;

    MarcherInfo& info;
    // grid spacing (for each axis)
    double const* const dx;
    // maximum order of the finite difference used
    int const order;

    // the shift to apply in every dimension to get a points neighbor (signed since we sometimes multiply this with -1)
    ssize* shift;
    // flags for each grid-point
    char* flags;

    // storage containers for solve_quadratic
    double *alpha_sq, *beta;
    bool* skip;
    auto solve_quadratic(usize const x, double* const tau) const -> double;

private:
    // inverse squared grid-spacing in each dimension (dx^2), used in solve_quadratic
    double* dx_sq_inv;

    void initialize(usize const x0, double* const tau);

public:
    Marcher(double const* const c, MarcherInfo& info, double const* dx, int const order);

    virtual ~Marcher();

    // virtual in the base class tells the compiler to do a late bind (see
    // http://www.willemer.de/informatik/cpp/cppvirt.htm)
    virtual void solve(usize const x0, double* const tau);
};
