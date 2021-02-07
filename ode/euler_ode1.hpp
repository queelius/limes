/**
 * First-order ODE solvers and the initial value problem (IVP).
 * 
 * The initial value problem is, given a derivative y' = f(t,y)
 * of an unknown function y(t) and an initial value y(t0) = y0,
 * find y(t).
 *
 * Note: if f(t,y) is a constant with respect to y, then y(t,y) = f(t),
 *       which can be found using normal integration, e.g., if the
 *       antiderivative F(t) := S f(t) dt can be found, then
 *           y(t) = F(t) + c
 *       where c is the constant of integration. To find the value of c,
 *       we use the IVP
 *           y(t0) = y0 => F(t0) + c = y0 => c = y0 - F(t0),
 *       thus
 *           y(t) = F(t) + y0 - F(t0).
 * 
 * If f(t,y) is a constant with respect to t, then f(t,y) = f(y).
 * 
 * We are given that y' = f(y), thus y = 
 * 
 * 
 * f(t) 
 */ 

#pragma once
#include <memory>
using std::forward;

namespace alex::math
{
    template <typename T>
    struct euler_ode1
    {
        using value_type = T;

        struct result_type
        {
            T operator()() const { return value; };
            T value;
            T error_estimate;
            int iterations;
        };

        T default_step_size() const { return 1e-3; };

        template <typename F>
        auto operator()(F&& f, T t0, T y0, T t) const
        {
            return operator()(forward<F>(f), t0, y0, t, default_step_size());
       	}

        template <typename F>
        auto operator()(F&& f, T t0, T y0, T t, T h) const
        {
            int i = 0;
            while (t0 < t - h)
            {
                y0 += f(t0, y0) * h;
                t0 += h;
                ++i;
            }
            if (t - t0 > 0)
            {
                ++i;
                y0 += f(t0, y0) * (t - t0);
            }
            return result_type{y0, h, i};
       	}
    };
}
