#include <vector>
#include <queue>
#include <math.h>
#include <stdexcept>
#include <numeric>
#include <functional>
#include "marcher.hpp"
#include "heap.cpp"


MarcherInfo::MarcherInfo(const int ndim, const usize *shape):
    ndim(ndim),
    shape(shape),
    size(std::accumulate(shape, shape+ndim, 1, std::multiplies<usize>()))
{}

void SensitivityInfo::store_sequence(const usize x)
{
    sequence[x] = counter++;
}

void SensitivityInfo::store_order(const int dim, const usize x, const char o)
{
    orders[dim*size + x] = o;
}


Marcher::Marcher(const double *const c, MarcherInfo &info, const double *dx, const int order):
    c(c),
    info(info),
    dx(dx),
    order(order)
{
    dx_sq_inv = new double[info.ndim];
    shift = new ssize[info.ndim];

    int cshift = 1;
    for (int i = info.ndim - 1; i >= 0; i--)
    {
        dx_sq_inv[i] = 1.0 / pow(dx[i], 2);

        shift[i] = cshift;
        cshift *= info.shape[i];
    }

    flags = new char[info.size];

    alpha_sq = new double[info.ndim];
    beta = new double[info.ndim];
    skip = new bool[info.ndim];
}

Marcher::~Marcher()
{
    delete[] shift;
    delete[] flags;

    delete[] dx_sq_inv;
    delete[] alpha_sq;
    delete[] beta;
    delete[] skip;
}

void Marcher::initialize(const usize x_s, double *tau)
{
    for (usize i = 0; i < info.size; i++)
    {
        flags[i] = UNKNOWN;
        tau[i] = INF;
    }

    tau[x_s] = 0;
}

void Marcher::solve(const usize x_s, double *const tau)
{
    initialize(x_s, tau);

    for (int d=0; d < info.ndim; d++)
        info.store_order(d, x_s, 0);

    auto heap_comp = [&tau](const usize e1, const usize e2){ return tau[e1] < tau[e2]; };
    Heap<decltype(heap_comp)> front(heap_comp, info.size);
    front.push(x_s);

    // list of points with minimal tau-values for each while-iteration
    std::vector<usize> minima;

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

        for (usize x_i : minima)
        {
            // store sequence number
            info.store_sequence(x_i);

            // now we check all neighbors of x_i
            usize rem = x_i;
            for (int d = 0; d < info.ndim; d++)
            {
                // index in the current axis
                usize dim_i = rem / shift[d];
                rem -= dim_i * shift[d];

                usize x_n = x_i - shift[d];
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
                if (dim_i < info.shape[d] - 1 && flags[x_n] != KNOWN)
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

double Marcher::solve_quadratic(const usize x, double *const tau) const
{
    usize rem = x;
    for (int d = 0; d < info.ndim; d++)
    {
        // index in the current axis
        usize dim_i = rem / shift[d];
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
        if (dim_i < info.shape[d] - 1 && flags[x + shift[d]] == KNOWN && tau[x + shift[d]] < tau_n)
        {
            tau_n = tau[x + shift[d]];
            direction = 1;
        }

        // only include summands that are valid, skip the rest
        if (direction != 0)
        {
            skip[d] = false;

            // valid known second neighbor + second order condition -> use second order
            if (order >= 2 && ((direction == -1 && dim_i > 1) || (direction == 1 && dim_i < info.shape[d] - 2)) && flags[x + 2 * direction * shift[d]] == KNOWN
                && tau_n > tau[x + 2 * direction * shift[d]])
            {
                alpha_sq[d] = nine_fourths * dx_sq_inv[d];
                beta[d] = one_third * (4.0 * tau_n - tau[x + 2 * direction * shift[d]]);
                info.store_order(d, x, 2*direction);
            }
            // otherwise fall back to first order
            else
            {
                alpha_sq[d] = dx_sq_inv[d];
                beta[d] = tau_n;
                info.store_order(d, x, direction);
            }
        }
        else
        {
            skip[d] = true;
            info.store_order(d, x, 0);
        }
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

        for (int d = 0; d < info.ndim; d++)
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

    // a is always positive, so the '+' solution is larger (causality)
    return (-b + sqrt(disc)) / (2.0 * a);
}
