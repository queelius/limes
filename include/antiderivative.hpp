
/**
 * An antidervative of f is integrate(f,a,x) + k for any a and k.
 * We are free to instantiate these values with any value.
 *
 * By default, we let a = 0 and k = 0.
 *
 * F models the concept of a univariate function to be integrated.
 * N models the concept of a numerical integrator.
 */
template <typename F, typename N>
struct antideriv_of
{
    F f;
    N nint;
    value_type a;
    value_type k;

    using value_type = N::value_type;

    antideriv_of(F f, N nint, double a = 0, double k = 0) :
        f(f), nint(nint), a(a), k(k) {}

    auto operator()(value_type x) const
    {
        return a < x ? k + nint(f,a,x)() : k - nint(f,x,a)();
    }
};

template <typename F, typename N>
auto deriv(antideriv_of<F,N> const & f)
{
    return f.f;
}

template <typename F, typename N = gauss_quad_integrator<kbn_sum<double>>>
auto antideriv(F const & f, N nint = N(), typename N::value_type a = 0.)
{
    return antideriv_of{f,nint,a}
}


/**
 * There are an infinite set of antiderivatives of a function f,
 *     F(a) := { integrate(f,x,a) + k | k in R },
 * where a is some number in the domain of f.
 *
 * If we wish to pick from F(a) the antiderivative that satisfies
 *     g(x0) = y0,
 * then we may do so by taking *any* antidervative of f in F(a) and adding
 * an appropriate constant (assuming that x0 is in domain of f and the
 * maximum-width interval that is a subset of the domain of f contains number a
 * and x0).
 *
 * We take the liberty of assuming that a typical value we apply the
 * antiderivative to is near x0, and so we reparameterize to F(x0) and solve
 * for the constant.
 */
template <typename F, typename N>
auto solve(antideriv_of<F,N> F, double x0, double y0)
{
    F.a = x0;
    F.k = y0;
    return F;
}

template <typename F, typename N>
auto solve_antideriv(F f, double x0, double y0, N nint = N())
{
    return antideriv_of<F,N>{ f, nint, x0, y0 }
}
