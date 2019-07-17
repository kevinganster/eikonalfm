#pragma once

#include <limits>

// memory leak detection
#ifdef _DEBUG
	#include "debugtests.h"
#endif


const double INF = std::numeric_limits<double>::infinity();

const char UNKNOWN = 0;
const char KNOWN = 1;
const char FRONT = 2;

// we don't want the function names to get 'mangled' by the compiler
extern "C"
{

class Marcher
{
protected:
	// slowness field
	const double* const c;

	// number of dimensions (of c)
	const int ndim;
	// shape (of c)
	const long* const shape;
	// total number of grid points
	size_t size;
	// grid spacing (for each axis)
	const double* const dx;
	// order of the finite difference used
	const int order;

	// the shift to apply in every dimension to get a points neighbor
	long* shift;

	// flags for each grid-point
	char* flags;

	virtual double solve_quadratic(const size_t x, double* const tau) const;

private:
	// inverse squared grid-spacing in each dimension (dx^2)
	double* dx_sq_inv;

	void initialize(const size_t x0, double* const tau);

public:
	Marcher(const double* const c, const int ndim, const long* const shape, const double* const dx, const int order);

	virtual ~Marcher();

	// virtual in the base class tells the compiler to do a late bind (see http://www.willemer.de/informatik/cpp/cppvirt.htm)
	virtual void solve(const size_t x0, double* const tau);
};

}