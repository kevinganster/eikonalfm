#include <vector>
#include <math.h>
#include <stdexcept>
#include "factoredmarcher.hpp"
#include "heap.cpp"


#define tau(x) (tau0[x] * tau1[x])

FactoredMarcher::FactoredMarcher(const double* const c, const int ndim, const size_t* const shape, const double* const dx, const int order) :
	Marcher(c, ndim, shape, dx, order)
{
}

const ptrdiff_t* FactoredMarcher::to_vector(size_t i)
{
    ptrdiff_t *x = new ptrdiff_t[ndim];
	size_t rem = i;

	for (int d = 0; d < ndim; d++)
	{
		x[d] = rem / shift[d];
		rem -= x[d] * shift[d];
	}

	return x;
}

void FactoredMarcher::initialize(double* const tau0, double* const tau1, const size_t x_s, const ptrdiff_t* const x_s_v)
{
	// iterating vector (see below)
    ptrdiff_t *x = new ptrdiff_t[ndim];
	for (int d = 0; d < ndim; d++)
		x[d] = 0;

	for (size_t i = 0; i < size; i++)
	{
		flags[i] = UNKNOWN;

		double distance = 0;
		for (int d = 0; d < ndim; d++)
		{
			distance += pow((double)((x[d] - x_s_v[d]) * dx[d]), 2);
		}

		tau0[i] = sqrt(distance);
		tau1[i] = INF;

		// increase vector-coordinates
		for (int d = ndim - 1; d >= 0; d--)
		{
			if (++x[d] >= (ptrdiff_t)shape[d])
				x[d] = 0;
			else
				break;
		}
	}

	delete[] x;
	tau1[x_s] = 1.0 / c[x_s];
}

#include <iostream>
void FactoredMarcher::solve(const size_t x_s, double* const tau1)
{
	double* tau0 = new double[size];
	const ptrdiff_t* const x_s_v = to_vector(x_s);
	initialize(tau0, tau1, x_s, x_s_v);

	auto heap_comp = [&tau0, &tau1](const size_t e1, const size_t e2) { return tau(e1) < tau(e2); };
	Heap<decltype(heap_comp)> front(heap_comp, size);
	front.push(x_s);

	// list of points with minimal tau-values for each while-iteration
	std::vector<size_t> minima;

	// main loop for fast marching
	while (!front.empty())
	{
		minima.clear();

		// add the point with minimal tau-value in front
		minima.push_back(front.top());
		// mark it as known
		flags[front.top()] = KNOWN;
		double value = tau(front.top());
		// remove it from the heap
		front.pop();

		// append all points that have the same (minimal) tau-value -> front.push will get faster by this
		while (true)
		{
			if (!front.empty() && tau(front.top()) == value)
			{
				minima.push_back(front.top());
				flags[front.top()] = KNOWN;
				front.pop();
			}
			else
				break;
		}

		for (size_t x_i : minima)
		{
			// now we check all neighbors of x_i
			size_t rem = x_i;
			for (int d = 0; d < ndim; d++)
			{
				// index in the current axis
				size_t dim_i = rem / shift[d];
				rem -= dim_i * shift[d];

				size_t x_n = x_i - shift[d];
				// valid neighbor to the 'left'
				if (dim_i > 0 && flags[x_n] != KNOWN)
				{
					tau1[x_n] = solve_quadratic(tau0, tau1, x_s_v, x_n);
					// add neighbor to front-heap if it is not yet on it
					if (flags[x_n] == UNKNOWN)
					{
						front.push(x_n);
						flags[x_n] = FRONT;
					}
					// otherwise update the heap
					else
						front.update(x_n);
				}

				x_n = x_i + shift[d];
				// valid neighbor to the 'right'
				if (dim_i < shape[d] - 1 && flags[x_n] != KNOWN)
				{
					tau1[x_n] = solve_quadratic(tau0, tau1, x_s_v, x_n);
					// add neighbor to front-heap if it is not yet on it
					if (flags[x_n] == UNKNOWN)
					{
						front.push(x_n);
						flags[x_n] = FRONT;
					}
					// otherwise update the heap
					else
						front.update(x_n);
				}
			}
		}
	}

	delete[] tau0;
	delete[] x_s_v;
}

double FactoredMarcher::solve_quadratic(const double* const tau0, const double* const tau1, const ptrdiff_t *const x_s, const size_t x)
{
	size_t rem = x;
	for (int d = 0; d < ndim; d++)
	{
		// index in the current axis
		size_t dim_i = rem / shift[d];
		rem -= dim_i * shift[d];

		double tau_n = INF;
		char direction = 0;
		// choose the direction considering the non-factored first order approximation (tau_i-1 < tau_i+1 -> backward direction)
		// valid known neighbor to the 'left'
		if (dim_i > 0 && flags[x - shift[d]] == KNOWN)
		{
			tau_n = tau0[x - shift[d]] * tau1[x - shift[d]];
			direction = -1;
		}
		// valid known neighbor to the 'right' (with smaller value)
		if (dim_i < shape[d] - 1 && flags[x + shift[d]] == KNOWN && tau(x + shift[d]) < tau_n)
		{
			tau_n = tau(x + shift[d]);
			direction = 1;
		}

		// only include summands that are valid, skip the rest
		if (direction != 0)
		{
			skip[d] = false;

			double alpha_d, beta_d;
			// valid known second neighbor + second order condition -> use second order
			if (order == 2 && ((direction == -1 && dim_i > 1) || (direction == 1 && dim_i < shape[d] - 2)) && flags[x + 2 * direction * shift[d]] == KNOWN
				&& tau_n > tau(x + 2 * direction * shift[d]))
			{
				alpha_d = 3.0 * tau0[x] / (2.0 * dx[d]) - direction * (ptrdiff_t)(dim_i - x_s[d]) * dx[d] / tau0[x];
				beta_d = tau0[x] * (4.0 * tau1[x + direction * shift[d]] - tau1[x + 2 * direction * shift[d]]) / (2 * dx[d] * alpha_d);
			}
			// otherwise fall back to first order
			else
			{
				alpha_d = tau0[x] / dx[d] - direction * (ptrdiff_t)(dim_i - x_s[d]) * dx[d] / tau0[x];
				beta_d = tau0[x] * tau1[x + direction * shift[d]] / (dx[d] * alpha_d);
			}

			alpha_sq[d] = pow(alpha_d, 2);
			beta[d] = beta_d;
		}
		else
			skip[d] = true;
	}


	double c_base = -1.0 / pow(this->c[x], 2);
	double a, b, c;
	double disc;

	// the currently biggest "beta" and it's dimension
	double biggest;
	int biggest_d;

	do
	{
		a = 0.0, b = 0.0;
		c = c_base;

		biggest = N_INF;
		biggest_d = -1;

		for (int d = 0; d < ndim; d++)
		{
			if (!skip[d])
			{
				a += alpha_sq[d];
				b -= 2.0 * alpha_sq[d] * beta[d];
				c += alpha_sq[d] * pow(beta[d], 2);

				if (beta[d] > biggest)
				{
					biggest_d = d;
					biggest = beta[d];
				}
			}
		}

		if (biggest_d == -1)
			throw std::runtime_error("Negative discriminant in solve_quadratic.");

		disc = pow(b, 2) - 4.0 * a * c;
		skip[biggest_d] = true;
	} while (disc < 0);

	// a is always positive, so the '+' solution is larger
	return (-b + sqrt(disc)) / (2.0 * a);
}
