#pragma once

#include <limits>
#include <stddef.h>

// memory leak detection
//#ifdef _DEBUG
//	#include "debugtests.h"
//#endif

// "sizes" used (signed and unsigned)
using ssize = long;
using usize = unsigned long;

const double INF = std::numeric_limits<double>::infinity();
const double N_INF = std::numeric_limits<double>::lowest();

const char UNKNOWN = 0;
const char KNOWN = 1;
const char FRONT = 2;


class MarcherInfo
{
public:
	// number of dimensions
	const int ndim;
	// shape (of the mesh)
    const usize *const shape;
    // total number of grid points
    usize size;

    MarcherInfo(const int ndim, const usize *shape);

	virtual ~MarcherInfo(){}

	// No-Op functions for sensitivities in the base class
	virtual void store_sequence(const usize x){}
	virtual void store_order(const int dim, const usize x, const char o){}
	virtual usize *get_sequence(){ return nullptr; }
	virtual char *get_order(){ return nullptr; }
};

class SensitivityInfo : public MarcherInfo
{
public:
	// execution sequence
	usize counter = 0;
	usize *sequence;
	// finite difference order
	char *orders;

	SensitivityInfo(const int ndim, const usize *shape): MarcherInfo{ndim, shape}
	{
		sequence = new usize[size];
		orders = new char[ndim*size];
	}

	~SensitivityInfo() override
	{
		// delete[] sequence;
		// delete[] orders;
	}

	void store_sequence(const usize x) override;
	void store_order(const int dim, const usize x, const char o) override;

	auto get_sequence() -> usize * override
	{
		return sequence;
	}

	auto get_order() -> char * override
	{
		return orders;
	}
};


class Marcher
{
protected:
	// slowness field
	const double *const c;

	MarcherInfo &info;
    // grid spacing (for each axis)
    const double *const dx;
	// maximum order of the finite difference used
	const int order;

	// the shift to apply in every dimension to get a points neighbor (signed since we sometimes multiply this with -1)
	ssize *shift;
	// flags for each grid-point
	char *flags;

	// storage containers for solve_quadratic
	double *alpha_sq, *beta;
	bool *skip;
	auto solve_quadratic(const usize x, double *const tau) const -> double;

private:
	// inverse squared grid-spacing in each dimension (dx^2), used in solve_quadratic
	double *dx_sq_inv;

	void initialize(const usize x0, double *const tau);

public:
	Marcher(const double *const c, MarcherInfo &info, const double *dx, const int order);

	virtual ~Marcher();

	// virtual in the base class tells the compiler to do a late bind (see http://www.willemer.de/informatik/cpp/cppvirt.htm)
	virtual void solve(const usize x0, double *const tau);
};
