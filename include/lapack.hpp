#ifndef TCM_LAPACK_HPP
#define TCM_LAPACK_HPP


#include <utils.hpp>
#include <hermitian.hpp>
#include <general.hpp>
#include <matrix.hpp>


namespace tcm {


namespace lapack {


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
template<class _Matrix, class _Vector>
inline
auto heev( _Matrix& A, _Vector& W
         , bool compute_eigenvectors ) -> void
{
	auto const N = A.height();

	assert(is_square(A));
	assert(is_column(W));
	assert(W.height() == N);

	lapack::heev( N
	            , A.data(), A.ldim()
	            , W.data()
	            , compute_eigenvectors );
}


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
template<class _Matrix, class _Vector>
auto heevr(_Matrix& A, _Vector& W) -> void
{
	auto const N = A.height();

	assert(is_square(A));
	assert(is_column(W));
	assert(W.height() == N);

	lapack::heevr( N
	             , A.data(), A.ldim()
	             , W.data()
	             , static_cast<typename _Matrix::value_type*>(nullptr), 1 );
}


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
template<class _Matrix1, class _Matrix2, class _Vector>
auto heevr(_Matrix1& A, _Vector& W, _Matrix2& Z) -> void
{
	auto const N = A.height();

	assert(is_square(A));
	assert(is_square(Z));
	assert(is_column(W));
	assert(W.height() == N);
	assert(Z.height() == N);

	lapack::heevr( N
	             , A.data(), A.ldim()
	             , W.data()
	             , Z.data(), Z.ldim() );
}


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
template<class _Matrix, class _Vector>
auto geev(_Matrix& A, _Vector& W) -> void
{
	auto const N = A.height();

	assert(is_square(A));
	assert(is_column(W));
	assert(W.height() == N);

	lapack::geev( N
	            , A.data(), A.ldim()
	            , W.data()
	            , static_cast<typename _Vector::value_type*>(nullptr), 1
	            , static_cast<typename _Vector::value_type*>(nullptr), 1
	            );
}


///////////////////////////////////////////////////////////////////////////////
///
///////////////////////////////////////////////////////////////////////////////
template<class _Matrix1, class _Matrix2, class _Vector>
auto geev(_Matrix1& A, _Vector& W, _Matrix2& Z) -> void
{
	auto const N = A.height();

	assert(is_square(A));
	assert(is_square(Z));
	assert(is_column(W));
	assert(W.height() == N);
	assert(Z.height() == N);

	lapack::geev( N
	            , A.data(), A.ldim()
	            , W.data()
	            , static_cast<typename _Vector::value_type*>(nullptr), 1
	            , Z.data(), Z.ldim()
	            );
}



} // namespace lapack

} // namespace tcm


#endif // TCM_LAPACK_HPP
