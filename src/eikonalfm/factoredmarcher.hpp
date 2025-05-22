#pragma once

#include "marcher.hpp"

class FactoredMarcher : public Marcher
{
private:
	const ssize* to_vector(usize i);
	void initialize(double* const tau0, double* const tau1, const usize x0, const ssize* const x0_v);

protected:
	double solve_quadratic(const double* const tau0, const double* const tau1, const ssize* const x0, const usize x);

public:
	FactoredMarcher(const double* const c, MarcherInfo& info, const double* const dx, const int order);

	void solve(const usize x0, double* const tau1);
};