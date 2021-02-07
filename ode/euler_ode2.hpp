/**
 * Second-order ODE initial value problem solver (IVP).
 * 
 * The initial value problem is, given
 *     y'' = f(t,y',y)
 * of an unknown function y(t) and an initial value
 *     y(t0)  = y0,
 *     y'(t0) = m0,
 * find y(t).
 * 
 * We transform this into a system of first-order equations.
 * Let z = y', then, z' = y''. Thus,
 *     y'    = f2(t,y)   = z
 *     z'    = f1(t,z,y) = f(t,z,y)
 *     y(t0) = y0,
 *     z(t0) = m0.
 * 
 * Each iteration of the numerical algorithm is of the form
 *     y{i+1} = yi + g(zi)
 *     z{i+1} = zi + g(f(ti,zi,yi))
 * where g is the step-algorithm. For instance, if using 
 * euler's method,
 *     y{i+1} = yi + zi . h
 *     z{i+1} = zi + f(ti,zi,yi) . h
 * Repeat this process k times such that t = k . h.
 */ 

#pragma once

#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

namespace alex::math
{
	template <typename T>
    struct euler_ode2
    {
		using value_type = T;

        struct result_type
        {
			T operator()() const { return y; };
            T y;
            T m;
            T error_estimate;
            int iterations;
        };

		T default_step_size() const { return 1e-3; };

		template <typename F>
        result_type ivp(F&& f, T t0, T y0, T t) const
        {
            return ivp(std::forward<F>(f), t0, y0, t, default_step_size());
       	};

        template <typename F>
        result_type bvp(F&& f, T t0, T y0, T t1, T y1, T lb, T ub, T h) const
        {
            // implemented via the shooting method using binary search
            T m;
            result_type s;
            int i = 0;
            while (std::abs(ub - lb) > 1e-12)
            {
                m = (lb + ub) / T(2);
                s = ivp(f, t0, y0, m, t1, h);

                i += s.iterations;
                if (s() < y1)
                    lb = m;
                else
                    ub = m;
            }
            return result_type{s(), m, 0, i};
        };

        template <typename F>
        result_type linear_bvp(F&& f, T t0, T y0, T t1, T y1, T lb, T ub, T h) const
        {
            // implemented via the interpolation method
            auto s = ivp(f, t0, y0, lb, t1, h);
            auto q0 = s.y;
            int i = s.iterations;

            s = ivp(f, t0, y0, ub, t1, h);
            auto q1 = s.y;
            i += s.iterations;
            
            auto m = lb + (ub - lb) / (q1 - q0) * (y1 - q0);

            s = ivp(f, std::move(t0), std::move(y0), m, std::move(t1), std::move(h));
            i += s.iterations;

            return result_type{std::move(s.y), std::move(m), std::move(s.error_estimate), i};
        };


		template <typename F>
        result_type ivp(F&& f, T t0, T y0, T m0, T t, T h) const
        {
            int i = 0;
            while (t0 < t - h)
            {
                auto y_old = y0;
                y0 += m0 * h;
                m0 += f(t0, std::move(y_old), m0) * h;
                
                t0 += h;
                ++i;
            }
            if (t - t0 > 0)
            {
                auto y_old = y0;
                y0 += m0 * (t - t0);
                m0 += f(t0, std::move(y_old), m0) * (t - t0);
                ++i;
            }
            return result_type{std::move(y0), std::move(m0), std::move(h), i};
       	};
    };
}
