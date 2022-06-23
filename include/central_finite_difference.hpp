#pragma once

#include <utility>

template <typename T>
struct differential_result
{
	T difference;
	T error_estimate;
	int evaluations;
};

// computes the numerical derivative of a function f at x
// using central finite difference coefficients where the
// accuracy is of order h^8.
//
// note that this computes up to 4h distance away from x, so if
// x is very near an asymptotic or undefined region, it may fail.
template <typename X, typename F>
auto central_finite_difference(
	F f,
	X x,
	X h = 1e-3)
{
	using Y = decltype(f(x));

	Y const A = .00357143;
	Y const B = .0380952;
	Y const C = .2;
	Y const D = .8;

	return differential_result<Y>{
		(A * f(x - X(4) * h) -
		 B * f(x - X(3) * h) +
		 C * f(x - X(2) * h) -
		 D * f(x - h) +
		 D * f(x + h) -
		 C * f(x + X(2) * h) +
		 B * f(x + X(3) * h) -
		 A * f(x + X(4) * h)) / static_cast<Y>(h),
		 h*h*h*h, 8};
}

template <typename X, typename F>
auto central_finite_difference_2nd(
	F f,
	X x,
	X h = 1e-3)
{
	using Y = decltype(f(x));

	Y const A = -0.00178571428;
	Y const B = 0.02539682539;
	Y const C = -0.2;
	Y const D = 1.6;
	Y const E = -2.84722222222;

	return differential_result<Y>
	{
		(A * f(x - X(4) * h) +
		 B * f(x - X(3) * h) +
		 C * f(x - X(2) * h) +
		 D * f(x - h) +
		 E * f(x) +
		 D * f(x + h) +
		 C * f(x + X(2) * h) +
		 B * f(x + X(3) * h) +
		 A * f(x + X(4) * h)) / static_cast<Y>(h * h),
		 h * h * h, 9
	};
}

template <typename X, typename F>
auto central_finite_difference_3rd(
	F f,
	X x,
	X h = 1e-3)
{
	using Y = decltype(f(x));

	Y const A = 0.02916666666;
	Y const B = 0.3;
	Y const C = 1.40833333333;
	Y const D = 2.03333333333;

	return differential_result<Y>
	{
		(-A * f(x - X(4) * h) +
		 B * f(x - X(3) * h) -
		 C * f(x - X(2) * h) +
		 D * f(x - h) -
		 D * f(x + h) +
		 C * f(x + X(2) * h) -
		 B * f(x + X(3) * h) +
		 A * f(x + X(4) * h)) / static_cast<Y>(h * h * h),
		 h * h * h, 8
	};
}

