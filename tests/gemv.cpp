#include <iostream>
#include <iomanip>
#include <cassert>
#include <unordered_map>

#include <blas.hpp>

using namespace tcm;


template<class _T>
auto apply_gemv(std::size_t const N, std::size_t const M) -> void
{
	Matrix<_T> A{N, M};
	Matrix<_T> V{M, 1};
	
	std::cin >> A >> V;

	Matrix<_T> Y{N, 1}; 
	blas::gemv( blas::Operator::None
	          , _T{1.0}, A, V
	          , _T{0.0}, Y );

	std::cout << std::setprecision(20) << Y << '\n';
}



int main(int argc, char** argv)
{
	std::unordered_map< std::string
	                  , void (*)(std::size_t const, std::size_t const) 
	                  > func_map;
	func_map["float"]          = &apply_gemv<float>;
	func_map["double"]         = &apply_gemv<double>;
	func_map["complex-float"]  = &apply_gemv<std::complex<float>>;
	func_map["complex-double"] = &apply_gemv<std::complex<double>>;

	assert(argc == 4);
	const auto N = static_cast<std::size_t>(std::stoi(argv[2]));
	const auto M = static_cast<std::size_t>(std::stoi(argv[3]));
	
	func_map.at(argv[1])(N, M);
	return 0;
}
