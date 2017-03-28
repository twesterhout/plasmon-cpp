#include <iostream>
#include <iomanip>
#include <cassert>
#include <map>

#define DO_MEASURE

#include <matrix.hpp>
#include <lapack.hpp>


using namespace tcm;


template<class T>
auto apply_geev(std::size_t const N) -> void
{
	Matrix<T, 64> A{N, N};
	Matrix<std::complex<utils::Base<T>>, 64> W{N, 1};
	
	std::cin >> A;

	lapack::geev(A, W);
	std::sort( W.begin_column(0), W.end_column(0)
	         , [](auto x, auto y) { return std::real(x) < std::real(x); } );
	std::cout << std::setprecision(20) << W << '\n';
}


int main(int argc, char** argv)
{
	std::map< std::string
	        , void (*)(std::size_t const) > func_map;

	// func_map["float"]          = &apply_geev<float>;
	// func_map["double"]         = &apply_geev<double>;
	func_map["complex-float"]  = &apply_geev<std::complex<float>>;
	func_map["complex-double"] = &apply_geev<std::complex<double>>;

	assert(argc == 3);
	const auto N = static_cast<std::size_t>(std::stoi(argv[2]));

	func_map.at(argv[1])(N);

	timing::report(std::cerr);
	return 0;
}

