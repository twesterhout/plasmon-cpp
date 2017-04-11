#include <iostream>
#include <iomanip>
#include <cassert>
#include <map>

#include <matrix.hpp>
#include <blas.hpp>

using namespace tcm;


template<class T>
auto apply_gemm( std::size_t const N
               , std::size_t const M
	       , std::size_t const K) -> void
{
	Matrix<T> A{N, M};
	Matrix<T> B{M, K};
	Matrix<T> C{N, K};
	
	std::cin >> A >> B;

	blas::gemm( blas::Operator::None, blas::Operator::None
	          , T{1.0}, A, B
		  , T{0.0}, C
		  );

	std::cout << std::setprecision(20) << C << '\n';
}


int main(int argc, char** argv)
{
	std::map< std::string
	        , void (*)( std::size_t const
	                  , std::size_t const
			  , std::size_t const) > func_map;

	func_map["float"]          = &apply_gemm<float>;
	func_map["double"]         = &apply_gemm<double>;
	func_map["complex-float"]  = &apply_gemm<std::complex<float>>;
	func_map["complex-double"] = &apply_gemm<std::complex<double>>;

	assert(argc == 5);
	const auto N = static_cast<std::size_t>(std::stoi(argv[2]));
	const auto M = static_cast<std::size_t>(std::stoi(argv[3]));
	const auto K = static_cast<std::size_t>(std::stoi(argv[4]));

	func_map.at(argv[1])(N, M, K);

	return 0;
}

