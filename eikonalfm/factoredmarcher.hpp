#pragma once

#include "marcher.hpp"

class FactoredMarcher : public Marcher
{
private:
	const ptrdiff_t* to_vector(size_t i);
	void initialize(double *const tau0, double *const tau1, const size_t x0, const ptrdiff_t *const x0_v);

protected:
	double solve_quadratic(const double *const tau0, const double *const tau1, const ptrdiff_t *const x0, const size_t x);

public:
	FactoredMarcher(const double *const c, const int ndim, const size_t *const shape, const double *const dx, const int order);

	void solve(const size_t x0, double *const tau1);
};