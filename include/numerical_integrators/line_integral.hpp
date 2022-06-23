#pragma once

#include <utility>
#include <vector>
#include <random>
#include <algorithm>

namespace numerical_integration
{
    template <typename F, typename X, typename Y, typename T = long double>
    auto adaptive_line_integral(
        F f, X x, Y y,
        T tb, T te, T dt, T const eps = T(0.001))
    {
        using std::min;
        using std::abs;

        static int N = 2;

        auto sum = T(0);
        auto dp = (te - tb) / T(N);
        while (tb < te)
        {
            auto end = min(te, tb + dp);
            auto sum0 = line_integral(f, x, y, tb, end, dt);
            auto sum1 = line_integral(f, x, y, tb, end, dt / T(N));

            if (abs(sum0 - sum1) > eps && dt > eps)
                sum1 = adaptive_line_integral(f, x, y, tb, tb + dp, dt / T(N));

            sum += sum1;
            tb += dp;
        }

        return sum;
    }

    template <typename F, typename X, typename Y, typename T = long double>
    auto line_integral(
        F f, X x, Y y,
        T tb, T te, T dt)
    {
        using std::sqrt;
        using std::pow;

        auto sum = T(0);

        while (tb < te)
        {
            auto x0 = x(tb);
            auto y0 = y(tb);

            tb += dt;
            if (tb > te)
                tb = te;

            auto y1 = y(tb);
            auto x1 = x(tb);

            sum += (f(x0, y0) + f(x1, y1)) *
                sqrt(pow(x1 - x0, T(2)) + pow(y1 - y0, T(2)));
        }

        return T(0.5) * sum;
    }

    template <typename X, typename Y, typename T = long double>
    auto path_integral(X x, Y y, T tb, T te, T dt)
    {
        return line_integral([] (T,T) { return T(1); }, x, y, tb, te, dt);
    }
};
