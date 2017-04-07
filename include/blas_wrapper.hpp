#ifndef TCM_BLAS_WRAPPER_HPP
#define TCM_BLAS_WRAPPER_HPP

#include <type_traits>
#include <complex>
#include <cassert>

#include <boost/numeric/conversion/cast.hpp>

#include <config.hpp>
#include <utils.hpp>
#include <iterator.hpp>


///////////////////////////////////////////////////////////////////////////////
/// \file blas_wrapper.hpp
/// \brief Defines template versions of some BLAS routines.
///////////////////////////////////////////////////////////////////////////////



#define REGISTER_BLAS(general_f, s_f, d_f, c_f, z_f)                           \
	namespace {                                                                \
		template<class... Args>                                                \
        __attribute__((always_inline))                                         \
		inline                                                                 \
		auto general_f##_impl( utils::Type2Type<float>                         \
		                     , Args&&... args) noexcept                        \
		{ return s_f(std::forward<Args>(args)...); }                           \
		                                                                       \
		template<class... Args>                                                \
        __attribute__((always_inline))                                         \
		inline                                                                 \
		auto general_f##_impl( utils::Type2Type<double>                        \
		                     , Args&&... args) noexcept                        \
		{ return d_f(std::forward<Args>(args)...); }                           \
		                                                                       \
		template<class... Args>                                                \
        __attribute__((always_inline))                                         \
		inline                                                                 \
		auto general_f##_impl( utils::Type2Type<std::complex<float>>           \
		                     , Args&&... args) noexcept                        \
		{ return c_f(std::forward<Args>(args)...); }                           \
		                                                                       \
		template<class... Args>                                                \
        __attribute__((always_inline))                                         \
		inline                                                                 \
		auto general_f##_impl( utils::Type2Type<std::complex<double>>          \
		                     , Args&&... args) noexcept                        \
		{ return z_f(std::forward<Args>(args)...); }                           \
	}                                                                          \
                                                                               \
	template<class T, class... Args>                                           \
	__attribute__((always_inline))                                             \
	inline                                                                     \
	auto general_f(Args&&... args) noexcept                                    \
	{ return general_f##_impl( utils::Type2Type<T>{}                           \
	                         , std::forward<Args>(args)...);                   \
	}   







namespace tcm {

namespace import {


enum class Operator : char {None = 'N', T = 'T', H = 'C'};



#ifdef USING_INTEL_MKL
	using blas_int = long long int;
#else
#	ifdef USING_ATLAS
		using blas_int = int;
#	else
#		error "Need BLAS"
#	endif
#endif



// ============================================================================
//                                DOT PRODUCT                                  
// ============================================================================


extern "C" {

float sdot_
( blas_int const* N
, float const* X, blas_int const* INCX
, float const* Y, blas_int const* INCY );

double ddot_
( blas_int const* N
, double const* X, blas_int const* INCX
, double const* Y, blas_int const* INCY );

std::complex<float> cdotc_
( blas_int const* N
, std::complex<float> const* X, blas_int const* INCX
, std::complex<float> const* Y, blas_int const* INCY );

std::complex<double> zdotc_
( blas_int const* N
, std::complex<double> const* X, blas_int const* inc_X
, std::complex<double> const* Y, blas_int const* inc_Y );

} // extern "C"

REGISTER_BLAS(dotc, sdot_, ddot_, cdotc_, zdotc_)


template< class _T
        , class = std::enable_if_t
                  <   std::is_same<_T, float>() 
                   or std::is_same<_T, double>()
                   or std::is_same<_T, std::complex<float>>()
                   or std::is_same<_T, std::complex<double>>()
                  >
        >
inline
auto dot( std::size_t const n
        , _T const* X, blas_int const INCX
        , _T const* Y, blas_int const INCY
        ) -> _T
{
	if(n == 0) return _T{0.0};

	assert(X != nullptr and INCX != 0);
	assert(Y != nullptr and INCY != 0);
	/*
	utils::assert_valid( const_blas_iterator<_T>{X           , INCX}
	                   , const_blas_iterator<_T>{X + n * INCX, INCX} );
	utils::assert_valid( const_blas_iterator<_T>{Y           , INCY}
	                   , const_blas_iterator<_T>{Y + n * INCY, INCY} );
	*/

	auto const N = boost::numeric_cast<blas_int>(n);	
	auto const result = dotc<_T>(&N, X, &INCX, Y, &INCY);

	utils::assert_valid(result);
	return result;
}



// ============================================================================
//                                   AXPY                                      
// ============================================================================

extern "C" {

void saxpy_
( blas_int const* N, float const* A
, float const* X, blas_int const* INCX
, float      * Y, blas_int const* INCY );

void daxpy_
( blas_int const* N, double const* A
, double const* X, blas_int const* INCX
, double      * Y, blas_int const* INCY );

void caxpy_
( blas_int const* N, std::complex<float> const* A
, std::complex<float> const* X, blas_int const* INCX
, std::complex<float>      * Y, blas_int const* INCY );

void zaxpy_
( blas_int const* N, std::complex<double> const* A
, std::complex<double> const* X, blas_int const* INCX
, std::complex<double>      * Y, blas_int const* INCY );

} // extern "C"

REGISTER_BLAS(axpy, saxpy_, daxpy_, caxpy_, zaxpy_)

template< class _T
        , class = std::enable_if_t
                  <   std::is_same<_T, float>() 
                   or std::is_same<_T, double>()
                   or std::is_same<_T, std::complex<float>>()
                   or std::is_same<_T, std::complex<double>>()
                  >
        >
inline
auto axpy( std::size_t const n, _T const a
         , _T const* X, blas_int const INCX
         , _T      * Y, blas_int const INCY ) -> void
{
	if(n == 0) return;

	assert(X != nullptr and INCX != 0);
	assert(Y != nullptr and INCY != 0);
	/*
	assert_valid( const_blas_iterator<_T>{X           , INCX}
	            , const_blas_iterator<_T>{X + n * INCX, INCX} );
	assert_valid( const_blas_iterator<_T>{Y           , INCY}
	            , const_blas_iterator<_T>{Y + n * INCY, INCY} );
	*/
	auto const N = boost::numeric_cast<blas_int>(n);	
	
	axpy<_T>(&N, &a, X, &INCX, Y, &INCY);

	/*
	assert_valid( const_blas_iterator<_T>{Y           , INCY}
	            , const_blas_iterator<_T>{Y + n * INCY, INCY} );
	*/
}






// ============================================================================
//                         MATRIX-VECTOR MULTIPLICATION                        
// ============================================================================

extern "C" {

void sgemv_
( char const* TRANSA
, blas_int const* M, blas_int const* N
, float const* ALPHA, float const* A, blas_int const* LDA
, float const* X, blas_int const* INCX
, float const* BETA, float* Y, blas_int const* INCY );

void dgemv_
( char const* TRANSA
, blas_int const* M, blas_int const* N
, double const* ALPHA, double const* A, blas_int const* LDA
, double const* X, blas_int const* INCX
, double const* BETA, double* Y, blas_int const* INCY );

void cgemv_
( char const* TRANSA
, blas_int const* M, blas_int const* N
, std::complex<float> const* ALPHA
, std::complex<float> const* A, blas_int const* LDA
, std::complex<float> const* X, blas_int const* INCX
, std::complex<float> const* BETA
, std::complex<float>* Y, blas_int const* INCY );

void zgemv_
( char const* TRANSA
, blas_int const* M, blas_int const* N
, std::complex<double> const* ALPHA
, std::complex<double> const* A, blas_int const* LDA
, std::complex<double> const* X, blas_int const* INCX
, std::complex<double> const* BETA
, std::complex<double>* Y, blas_int const* INCY );

} // extern "C"

REGISTER_BLAS(gemv, sgemv_, dgemv_, cgemv_, zgemv_)


template< class _T
        , class = std::enable_if_t
                  <   std::is_same<_T, float>() 
                   or std::is_same<_T, double>()
                   or std::is_same<_T, std::complex<float>>()
                   or std::is_same<_T, std::complex<double>>()
                  >
        >
auto gemv( Operator const op_A
         , std::size_t const m, std::size_t const n
         , _T const ALPHA, _T const* A, blas_int const LDA
         , _T const* X, blas_int const INCX
         , _T const BETA, _T* Y, blas_int const INCY ) -> void
{
	if(n == 0 or m == 0) return;

	assert(A != nullptr and LDA  != 0);
	assert(X != nullptr and INCX != 0);
	assert(Y != nullptr and INCY != 0);
	/*
	utils::assert_valid(A, A + LDA * n);
	utils::assert_valid
		( const_blas_iterator<_T>{X, INCX}
		, const_blas_iterator<_T>{ X + (op_A == Operator::None ? n : m) * INCX
		                         , INCX }
		);
	utils::assert_valid
		( const_blas_iterator<_T>{Y, INCY}
		, const_blas_iterator<_T>{ Y + (op_A == Operator::None ? m : n) * INCY
		                         , INCY }
		);
	*/

	auto const TRANS = static_cast<char>(op_A);
	auto const M     = boost::numeric_cast<blas_int>(m);
	auto const N     = boost::numeric_cast<blas_int>(n);

	assert(LDA >= M);

	gemv<_T>( &TRANS, &M, &N
	        , &ALPHA, A, &LDA, X, &INCX
	        , &BETA, Y, &INCY );

	/*
	utils::assert_valid
		( const_blas_iterator<_T>{Y, INCY}
		, const_blas_iterator<_T>{ Y + (op_A == Operator::None ? m : n) * INCY
		                         , INCY }
		);
	*/
}




// ============================================================================
//                           MATRIX-MATRIX MULTIPLICATION                      
// ============================================================================

extern "C" {

void sgemm_
( char const* TRANSA, char const* TRANSB
, blas_int const* M, blas_int const* N, blas_int const* K
, float const* ALPHA
, float const* A, blas_int const* LDA
, float const* B, blas_int const* LDB
, float const* BETA
, float* C, blas_int const* LDC );

void dgemm_
( char const* TRANSA, char const* TRANSB
, blas_int const* M, blas_int const* N, blas_int const* K
, double const* ALPHA
, double const* A, blas_int const* LDA
, double const* B, blas_int const* LDB
, double const* BETA
, double* C, blas_int const* LDC );

void cgemm_
( char const* TRANSA, char const* TRANSB
, blas_int const* M, blas_int const* N, blas_int const* K
, std::complex<float> const* ALPHA
, std::complex<float> const* A, blas_int const* LDA
, std::complex<float> const* B, blas_int const* LDB
, std::complex<float> const* BETA
, std::complex<float>* C, blas_int const* LDC );

void zgemm_
( char const* TRANSA, char const* TRANSB
, blas_int const* M, blas_int const* N, blas_int const* K
, std::complex<double> const* ALPHA
, std::complex<double> const* A, blas_int const* LDA
, std::complex<double> const* B, blas_int const* LDB
, std::complex<double> const* BETA
, std::complex<double>* C, blas_int const* LDC );

} // extern "C"

REGISTER_BLAS(gemm, sgemm_, dgemm_, cgemm_, zgemm_)

template< class _T
        , class = std::enable_if_t
                  <   std::is_same<_T, float>() 
                   or std::is_same<_T, double>()
                   or std::is_same<_T, std::complex<float>>()
                   or std::is_same<_T, std::complex<double>>()
                  >
        >
inline
auto gemm( Operator const op_A, Operator const op_B
         , std::size_t const m, std::size_t const n, std::size_t const k
         , _T const ALPHA, _T const* A, blas_int const LDA
         , _T const* B, blas_int const LDB
         , _T const BETA, _T* C, blas_int const LDC ) -> void
{
	if(m == 0 or n == 0 or k == 0) return;

	assert(A != nullptr and LDA != 0);
	assert(B != nullptr and LDB != 0);
	assert(C != nullptr and LDC != 0);
	/*
	utils::assert_valid(A, A + LDA * (op_A == Operator::None ? k : m));
	utils::assert_valid(B, B + LDB * (op_B == Operator::None ? n : k));
	utils::assert_valid(C, C + LDC * n);
	*/

	auto const M      = boost::numeric_cast<blas_int>(m);
	auto const N      = boost::numeric_cast<blas_int>(n);
	auto const K      = boost::numeric_cast<blas_int>(k);
	auto const TRANSA = static_cast<char>(op_A);
	auto const TRANSB = static_cast<char>(op_B);

	assert( LDA >= (op_A == Operator::None ? M : K) );
	assert( LDB >= (op_B == Operator::None ? K : N) );
	assert( LDC >= M );

	gemm<_T>( &TRANSA, &TRANSB, &M, &N, &K
	        , &ALPHA, A, &LDA, B, &LDB
	        , &BETA, C, &LDC );

	/*
	utils::assert_valid(C, C + LDC * n);
	*/
}



} // namespace import

} // namespace tcm








#endif // TCM_BLAS_WRAPPER_HPP
