//#define DEBUG
#include <vector>
#include <iostream>
#include "euler_ode1.hpp"
#include "euler_ode2.hpp"
//#include "rk4_ode1.hpp"
#include "rk4_ode2.hpp"
#include "../integration/kbn_accumulate.hpp"
#include <cmath>

int main()
{
    using acc = alex::math::kbn_accumulate<long double>;
    
    alex::math::rk4_ode2<acc> rk4;
    alex::math::euler_ode2<long double> euler2;
 
    // x'' + d x' + a x + b x^3 = r cos(w t)
    // x'' = -d x' - a x - b x^3 + r cos(w t)
    // y(0) = 1
    // y'(0) = 0
    //auto f = [d = .3, a = -1, b = 1, w = 1.2, r = .5](long double t, double y, double m)
    //{
    //    return -d*m - a*m - b*y*y*y + r * std::cos(w * t);
    //};
    // F = -k x - a v
    std::vector<bool> s(10,true);
    long double pulse = 1;
    auto f = [&s, &pulse, k = 2, a = 1](long double t, double y, double m)
    {
        auto p = -k * y - a * m;

        auto t_int = int(t);
        if (t > pulse && t_int < s.size() && s[t_int-1])
        {
            std::cout << "t = " << t << std::endl;
            std::cout << "t_int = " << t_int << std::endl;
            std::cout << "p = " << p << std::endl;
            pulse++;
            s[t_int-1] = false;
            p += 1e20;
        }

        return p;
    };
    //auto f = [](long double t, double y, double m) { return -9.8 - std::sqrt(2)*m + t*std::sin(3*t); };
    //long double t0 = 0;
    //long double y0 = 100;
    //long double y1 = 0;
    //long double t1 = 5;

    auto r1 = rk4.ivp(f, 0, 0, 10, 20, 1e-3);
    std::cout << r1.m << std::endl;
    //auto r3 = euler2.ivp(f, t0, y0, m0, t1, 1e-3);

    //std::cout.precision(20);
    //std::cout << "y(" << t1 << ") = " << r1.y << std::endl;
    //std::cout << "y'(0) = " << r1.m << std::endl;
    //std::cout << "---\n";
    //std::cout << "y(" << t1 << ") = " << r3.y << std::endl;
    //std::cout << "y'(0) = " << r3.m << std::endl;

    return EXIT_SUCCESS;
}