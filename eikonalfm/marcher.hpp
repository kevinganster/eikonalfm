#pragma once

#include <limits>
#include <stddef.h>

// memory leak detection
//#ifdef _DEBUG
//	#include "debugtests.h"
//#endif


const double INF = std::numeric_limits<double>::infinity();
const double N_INF = std::numeric_limits<double>::lowest();

const char UNKNOWN = 0;
const char KNOWN = 1;
const char FRONT = 2;


class Marcher
{
protected:
	// slowness field
	const double* const c;

	// number of dimensions (of c)
	const int ndim;
	// shape (of c)
	const size_t* const shape;
	// total number of grid points
	size_t size;
	// grid spacing (for each axis)
	const double* const dx;
	// order of the finite difference used
	const int order;

	// the shift to apply in every dimension to get a points neighbor (ptrdiff_t since we sometimes multiply this with -1)
	ptrdiff_t* shift;

	// flags for each grid-point
	char* flags;

	double solve_quadratic(const size_t x, double* const tau) const;

	// storage containers for solve_quadratic
	double* alpha_sq,* beta;
	bool* skip;

private:
	// inverse squared grid-spacing in each dimension (dx^2), used in solve_quadratic
	double* dx_sq_inv;

	void initialize(const size_t x0, double* const tau);

public:
	Marcher(const double* const c, const int ndim, const size_t* const shape, const double* const dx, const int order);

	virtual ~Marcher();

	// virtual in the base class tells the compiler to do a late bind (see http://www.willemer.de/informatik/cpp/cppvirt.htm)
	virtual void solve(const size_t x0, double* const tau);
};
