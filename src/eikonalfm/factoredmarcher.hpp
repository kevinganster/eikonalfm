#pragma once

#include "marcher.hpp"
#include <vector>

class FactoredMarcher: public Marcher {
private:
    auto to_vector(usize i) -> std::vector<usize>;
    void initialize(double* const tau0, double* const tau1, const usize x_s, std::vector<usize> const& x_s_v);

protected:
    auto
    solve_quadratic(double const* const tau0, double const* const tau1, std::vector<usize> const& x_s_v, usize const x)
        -> double;

public:
    FactoredMarcher(double const* const c, MarcherInfo& info, double const* const dx, int const order);

    void solve(usize const x_s, double* const tau1) override;
};
