#pragma once

#include <vector>
#include <limits>

namespace calckit
{
    struct five_point_stencil_method
    {
        template <typename F, typename T>
        auto operator()(F f, T t, T h)
        {
            return (T(.00357143)*f(t - T(4)*h) - T(.0380952)*f(t - T(3)*h) + T(.2)*f(t - T(2)*h) - T(.8)*f(t - h) +
                    T(.8)*f(t + h) - T(.2)*f(t + T(2)*h) + T(.0380952)*f(t + T(3)*h) - T(.00357143)*f(t + T(4)*h)) / h;
        };
    };
    
    // f : R^m -> R, x in R^m, h in R
    template <typename F, typename X, typename T,
              typename D = five_point_stencil_method>
    auto grad(F f, X x, T h)
    {
        auto g = x;
        int i = 0;
        for (auto v : x)
        {
            g[i] = D{}([f, x, i](auto t) mutable {
                x[i] = t;
                return f(x);
            }, v, h);
            ++i;
        }
        return g;
    }
}
