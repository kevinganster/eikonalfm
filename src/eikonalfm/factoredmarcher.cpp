#include "factoredmarcher.hpp"
#include "heap.cpp"
#include <cmath>
#include <stdexcept>


#define tau(x) (tau0[x] * tau1[x])

FactoredMarcher::FactoredMarcher(double const* const c, MarcherInfo& info, double const* const dx, int const order)
    : Marcher(c, info, dx, order) {}

auto FactoredMarcher::to_vector(usize i) -> std::vector<usize> {
    std::vector<usize> x(info.ndim);
    usize rem = i;

    for (int d = 0; d < info.ndim; d++) {
        x[d] = rem / shift[d];
        rem -= x[d] * shift[d];
    }

    return x;
}

void FactoredMarcher::initialize(
    double* const tau0, double* const tau1, usize const x_s, std::vector<usize> const& x_s_v
) {
    // iterating vector (see below)
    std::vector<usize> x(info.ndim, 0);

    for (usize i = 0; i < info.size; ++i) {
        flags[i] = UNKNOWN;

        double distance = 0;
        for (int d = 0; d < info.ndim; d++)
            distance += std::pow((static_cast<double>(x[d]) - x_s_v[d]) * dx[d], 2);

        tau0[i] = std::sqrt(distance);
        tau1[i] = INF;

        // increase vector-coordinates
        for (int d = info.ndim - 1; d >= 0; --d) {
            if (++x[d] >= info.shape[d])
                x[d] = 0;
            else
                break;
        }
    }

    tau1[x_s] = 1.0 / c[x_s];
}

void FactoredMarcher::solve(usize const x_s, double* const tau1) {
    double* tau0 = new double[info.size];
    auto const x_s_v = to_vector(x_s);
    initialize(tau0, tau1, x_s, x_s_v);

    for (int d = 0; d < info.ndim; ++d) {
        info.store_order(d, x_s, 0);
    }

    auto heap_comp = [&tau0, &tau1](const usize e1, const usize e2) { return tau(e1) < tau(e2); };
    Heap<decltype(heap_comp)> front(heap_comp, info.size);
    front.push(x_s);

    // list of points with minimal tau-values for each while-iteration
    std::vector<usize> minima;

    // main loop for fast marching
    while (!front.empty()) {
        minima.clear();

        // add the point with minimal tau-value in front
        minima.push_back(front.top());
        // mark it as known
        flags[front.top()] = KNOWN;
        double value = tau(front.top());
        // remove it from the heap
        front.pop();

        // append all points that have the same (minimal) tau-value -> front.push will get faster by this
        while (true) {
            if (!front.empty() && tau(front.top()) == value) {
                minima.push_back(front.top());
                flags[front.top()] = KNOWN;
                front.pop();
            } else
                break;
        }

        for (usize x_i : minima) {
            // store sequence number
            info.store_sequence(x_i);

            // now we check all neighbors of x_i
            usize rem = x_i;
            for (int d = 0; d < info.ndim; d++) {
                // index in the current axis
                usize dim_i = rem / shift[d];
                rem -= dim_i * shift[d];

                usize x_n = x_i - shift[d];
                // valid neighbor to the 'left'
                if (dim_i > 0 && flags[x_n] != KNOWN) {
                    tau1[x_n] = solve_quadratic(tau0, tau1, x_s_v, x_n);
                    // add neighbor to front-heap if it is not yet on it
                    if (flags[x_n] == UNKNOWN) {
                        front.push(x_n);
                        flags[x_n] = FRONT;
                    }
                    // otherwise update the heap
                    else
                        front.update(x_n);
                }

                x_n = x_i + shift[d];
                // valid neighbor to the 'right'
                if (dim_i < info.shape[d] - 1 && flags[x_n] != KNOWN) {
                    tau1[x_n] = solve_quadratic(tau0, tau1, x_s_v, x_n);
                    // add neighbor to front-heap if it is not yet on it
                    if (flags[x_n] == UNKNOWN) {
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
}

auto FactoredMarcher::solve_quadratic(
    double const* const tau0, double const* const tau1, std::vector<usize> const& x_s, usize const x
) -> double {
    usize rem = x;
    for (int d = 0; d < info.ndim; ++d) {
        // index in the current axis
        usize dim_i = rem / shift[d];
        rem -= dim_i * shift[d];

        double tau_n = INF;
        signed char direction = 0;
        // choose the direction considering the non-factored first order approximation (tau_i-1 < tau_i+1 -> backward
        // direction) valid known neighbor to the 'left'
        if (dim_i > 0 && flags[x - shift[d]] == KNOWN) {
            tau_n = tau0[x - shift[d]] * tau1[x - shift[d]];
            direction = -1;
        }
        // valid known neighbor to the 'right' (with smaller value)
        if (dim_i < info.shape[d] - 1 && flags[x + shift[d]] == KNOWN && tau(x + shift[d]) < tau_n) {
            tau_n = tau(x + shift[d]);
            direction = 1;
        }

        // only include summands that are valid, skip the rest
        if (direction != 0) {
            skip[d] = false;
            double alpha_d, beta_d;

            // valid known second neighbor + second order condition -> use second order
            if (order == 2 && ((direction == -1 && dim_i > 1) || (direction == 1 && dim_i < info.shape[d] - 2)) &&
                flags[x + 2 * direction * shift[d]] == KNOWN && tau_n > tau(x + 2 * direction * shift[d])) {
                // careful: long - ulong = ulong - long = ulong
                alpha_d =
                    3.0 * tau0[x] / (2.0 * dx[d]) - direction * (static_cast<double>(dim_i) - x_s[d]) * dx[d] / tau0[x];
                beta_d = tau0[x] * (4.0 * tau1[x + direction * shift[d]] - tau1[x + 2 * direction * shift[d]]) /
                         (2 * dx[d] * alpha_d);
                info.store_order(d, x, 2 * direction);
            }
            // otherwise fall back to first order
            else {
                alpha_d = tau0[x] / dx[d] - direction * (static_cast<double>(dim_i) - x_s[d]) * dx[d] / tau0[x];
                beta_d = tau0[x] * tau1[x + direction * shift[d]] / (dx[d] * alpha_d);
                info.store_order(d, x, direction);
            }

            alpha_sq[d] = alpha_d * alpha_d;
            beta[d] = beta_d;
        } else {
            skip[d] = true;
            info.store_order(d, x, 0);
        }
    }


    double c_base = -1.0 / std::pow(this->c[x], 2);
    double a, b, c;
    double disc;

    // the currently biggest "beta" and it's dimension
    double biggest;
    int biggest_d;

    do {
        a = 0.0, b = 0.0;
        c = c_base;

        biggest = N_INF;
        biggest_d = -1;

        for (int d = 0; d < info.ndim; d++) {
            if (!skip[d]) {
                a += alpha_sq[d];
                b -= 2.0 * alpha_sq[d] * beta[d];
                c += alpha_sq[d] * beta[d] * beta[d];

                if (beta[d] > biggest) {
                    biggest_d = d;
                    biggest = beta[d];
                }
            }
        }

        if (biggest_d == -1)
            throw std::runtime_error("Negative discriminant in solve_quadratic.");

        disc = b * b - 4.0 * a * c;
        skip[biggest_d] = true;
    } while (disc < 0);

    // a is always positive, so the '+' solution is larger
    return 0.5 * (-b + std::sqrt(disc)) / a;
}
