integration_test: univariate_integrator.hpp double_exponential_integrator.hpp integration_test.cpp gauss_quad_integrator.hpp kbn_accumulate.hpp numerical_univariate_integrator.hpp
	g++ -O3 -std=c++2a -Wall -o integration_test integration_test.cpp -I.
#-Wextra