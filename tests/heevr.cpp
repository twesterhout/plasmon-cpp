#include <iostream>
#include <iomanip>
#include <cassert>
#include <map>

#define DO_MEASURE

#include <matrix.hpp>
#include <lapack.hpp>

using namespace tcm;


template<class T>
auto apply_heevr(std::size_t const N) -> void
{
	Matrix<T> A{N, N};
	Matrix<utils::Base<T>> W{N, 1};
	Matrix<T> Z{N, N};
	
	std::cin >> A;

	lapack::heevr(A, W, Z);

	std::cout << std::setprecision(20) << W << '\n';
}


int main(int argc, char** argv)
{
	std::map< std::string
	        , void (*)(std::size_t const) > func_map;

	func_map["float"]          = &apply_heevr<float>;
	func_map["double"]         = &apply_heevr<double>;
	func_map["complex-float"]  = &apply_heevr<std::complex<float>>;
	func_map["complex-double"] = &apply_heevr<std::complex<double>>;

	assert(argc == 3);
	const auto N = static_cast<std::size_t>(std::stoi(argv[2]));

	func_map.at(argv[1])(N);

	timing::report(std::cerr);
	return 0;
}


