#pragma once

#include "marcher.hpp"

class FactoredMarcher : public Marcher
{
private:
	const long *to_vector(unsigned long i);
	void initialize(double *const tau0, double *const tau1, const unsigned long x0, const long *const x0_v);

protected:
	double solve_quadratic(const double *const tau0, const double *const tau1, const long *const x0, const unsigned long x);

public:
	FactoredMarcher(const double *const c, const int ndim, const unsigned long *const shape, const double *const dx, const int order);

	void solve(const unsigned long x0, double *const tau1);
};