#include <iostream>
#include <iomanip>
#include <cassert>
#include <unordered_map>

#include <matrix.hpp>
#include <blas.hpp>


template<class _T>
auto apply_dot(std::size_t const N) -> void
{
	tcm::Matrix<_T> X{N, 1};
	tcm::Matrix<_T> Y{N, 1};
	
	std::cin >> X;
	std::cin >> Y;

	std::cout << std::setprecision(20)
	          << tcm::blas::dot(X, Y) << '\n';
}


int main(int argc, char** argv)
{
	std::unordered_map< std::string
	                  , void (*)(std::size_t const) > func_map;
	func_map["float"]          = &apply_dot<float>;
	func_map["double"]         = &apply_dot<double>;
	func_map["complex-float"]  = &apply_dot<std::complex<float>>;
	func_map["complex-double"] = &apply_dot<std::complex<double>>;

	assert(argc == 3);
	const auto N = static_cast<std::size_t>(std::stoi(argv[2]));

	func_map.at(argv[1])(N);
	return 0;
}

