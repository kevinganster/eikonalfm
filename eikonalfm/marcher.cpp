#include <vector>
#include <queue>
#include <math.h>
#include "marcher.h"
#include "heap.cpp"

// we don't want the function names to get 'mangled' by the compiler
extern "C"
{

	Marcher::Marcher(const double* const c, const int ndim, const long* const shape, const double* const dx, const int order):
		c(c),
		ndim(ndim),
		shape(shape),
		dx(dx),
		order(order)
	{
		size = 1;
		dx_sq_inv = new double[ndim];
		shift = new long[ndim];

		for (int i = ndim - 1; i >= 0; i--)
		{
			dx_sq_inv[i] = 1.0 / pow(dx[i], 2);

			shift[i] = size;
			size *= shape[i];
		}

		flags = new char[size];
	}

	Marcher::~Marcher()
	{
		delete[] dx_sq_inv;
		delete[] shift;
		delete[] flags;
	}

	void Marcher::initialize(const size_t x0, double* tau)
	{
		for (long i = 0; i < size; i++)
		{
			flags[i] = UNKNOWN;
			tau[i] = INF;
		}

		tau[x0] = 0;
	}

	void Marcher::solve(const size_t x0, double* const tau)
	{
		initialize(x0, tau);

		auto heap_comp = [&tau](const size_t e1, const size_t e2){ return tau[e1] < tau[e2]; };
		Heap<decltype(heap_comp)> front(heap_comp, size);
		front.push(x0);

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
			double value = tau[front.top()];
			// remove it from the heap
			front.pop();

			// append all points that have the same (minimal) tau-value -> front.push will get faster by this
			bool done = false;
			while (!done)
			{
				if (!front.empty() && tau[front.top()] == value)
				{
					minima.push_back(front.top());
					flags[front.top()] = KNOWN;
					front.pop();
				}
				else
					done = true;
			}

			for (size_t x_i : minima)
			{
				// now we check all neighbors of x_i
				size_t rem = x_i;
				for (int d = 0; d < ndim; d++)
				{
					// index in the current axis
					long dim_i = rem / shift[d];
					rem -= dim_i * shift[d];

					size_t x_n = x_i - shift[d];
					// valid neighbor to the 'left'
					if (dim_i > 0 && flags[x_n] != KNOWN)
					{
						tau[x_n] = solve_quadratic(x_n, tau);
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
						tau[x_n] = solve_quadratic(x_n, tau);
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
	}

	const double nine_fourths = 9.0 / 4.0;
	const double one_third = 1.0 / 3.0;

	double Marcher::solve_quadratic(const size_t x, double* const tau) const
	{
		double a = 0.0, b = 0.0;
		double c = -1 / pow(this->c[x], 2);

		size_t rem = x;
		for (int d = 0; d < ndim; d++)
		{
			// index in the current axis
			long dim_i = rem / shift[d];
			rem -= dim_i * shift[d];

			double tau_n = INF;
			char direction = 0;
			// choose the direction considering the non-factored first order approximation (tau_i-1 < tau_i+1 -> backward direction)
			// valid known neighbor to the 'left'
			if (dim_i > 0 && flags[x - shift[d]] == KNOWN)
			{
				tau_n = tau[x - shift[d]];
				direction = -1;
			}
			// valid known neighbor to the 'right' (with smaller value)
			if (dim_i < shape[d] - 1 && flags[x + shift[d]] == KNOWN && tau[x + shift[d]] < tau_n)
			{
				tau_n = tau[x + shift[d]];
				direction = 1;
			}

			// only include summands that are valid, skip the rest
			if (direction != 0)
			{
				double alpha_d_sq, beta_d;
				// valid known second neighbor + second order condition -> use second order
				if (order == 2 && ((direction == -1 && dim_i > 1) || (direction == 1 && dim_i < shape[d] - 2)) && flags[x + 2 * direction * shift[d]] == KNOWN
					&& tau_n > tau[x + 2 * direction * shift[d]])
				{
					alpha_d_sq = nine_fourths;
					beta_d = one_third * (4.0 * tau_n - tau[x + 2 * direction * shift[d]]);
				}
				// otherwise fall back to first order
				else
				{
					alpha_d_sq = 1.0;
					beta_d = tau_n;
				}

				a += dx_sq_inv[d] * alpha_d_sq;
				b -= 2.0 * dx_sq_inv[d] * alpha_d_sq * beta_d;
				c += dx_sq_inv[d] * alpha_d_sq * pow(beta_d, 2);
			}
		}

		double disc = pow(b, 2) - 4.0 * a * c;

		if (disc < 0)
			throw std::runtime_error("Negative discriminant in solve_quadratic.");

		// a is always positive, so the '+' solution is larger
		return (-b + sqrt(disc)) / (2.0 * a);
	}
}