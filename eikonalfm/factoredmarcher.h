#pragma once

#include "marcher.h"

class FactoredMarcher : public Marcher
{
private:
	const long* to_vector(size_t i);
	void initialize(double* const tau0, double* const tau1, const size_t x0, const long* const x0_v, double* const tau);

protected:
	double solve_quadratic(const double* const tau0, const double* const tau1, const long* const x0, const size_t x, const double* const tau);

public:
	FactoredMarcher(const double* const c, const int ndim, const long* const shape, const double* const dx, const int order);

	void solve(const size_t x0, double* const tau);
};