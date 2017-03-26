#ifndef TCM_BLAS_HPP
#define TCM_BLAS_HPP


#include <matrix.hpp>
#include <blas_wrapper.hpp>
#include <benchmark.hpp>


///////////////////////////////////////////////////////////////////////////////
/// \file blas.hpp
/// \brief BLAS-like functionality
///////////////////////////////////////////////////////////////////////////////



namespace tcm {


namespace blas {

///////////////////////////////////////////////////////////////////////////////
/// \brief Enum class to denote BLAS matrix operator. 

/// It can be 
/// * `None`, i.e. \f$ A \mapsto A \f$;
/// * `T`, i.e. \f$ A \mapsto A^\text{T} \f$;
/// * `H`, i.e. \f$ A \mapsto A^\text{H} = \left(A^\text{T}\right)^* \f$.
///////////////////////////////////////////////////////////////////////////////
using Operator = import::Operator;

///////////////////////////////////////////////////////////////////////////////
/// \brief Typedef for integers used within BLAS. Allows for a more
/// <em>uniform</em> interface, as ATLAS may use `int`'s and Intel MKL
/// `long int`'s...
///////////////////////////////////////////////////////////////////////////////
using blas_int = import::blas_int;



// ============================================================================
// ||                                                                        ||
// ||                               DOT PRODUCT                              ||
// ||                                                                        ||
// ============================================================================


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculates dot product of two vectors.

/// This function is an abstraction from the ?DOT and ?DOTC functions in BLAS 
/// level 1. It dispatches to SDOT, DDOT, CDOTC or ZDOTC depending on the 
/// type of the matrices. 
/// Dot product is defined as
/// \f[ \langle X, Y \rangle := \sum_n X_n^* Y_n. \f]
/// 
/// \tparam _Vector1   Matrix/Vector type representing \f$\mathcal{F}^N\f$.
/// \tparam _Vector2   Matrix/Vector type representing \f$\mathcal{F}^N\f$.
///
/// Usually, _Vector1 and _Vector2 will just be %tcm::Matrix, but it's not
/// required. The following assumptions are made about _Vector1 and _Vector2:
/// 1) It has a value_type typedef, that is
///    `typename _Vector?::value_type` represents a type.
/// 2) is_row() and is_column() are defined for this type.
/// 3) `_Vector?::width()` returns the width of the vector if seen as a
///    matrix.
/// 4) `_Vector?::height()` returns the height of the vector if seen as a
///    matrix.
/// 5) `_Vector?::ldim()` returns the leading dimension of the vector.
/// 6) `_Vector?::data()` returns a pointer to the underlying array of
///    elements.
///
/// \param[in] X       First vector.
/// \param[in] Y       Second vector.
/// \return Dot product of X and Y.
/// \exception Nothrow guarantee.
///////////////////////////////////////////////////////////////////////////////
template<class _Vector1, class _Vector2>
inline
auto dot( _Vector1 const& X
        , _Vector2 const& Y ) noexcept -> typename _Vector1::value_type
{
	MEASURE;
	static_assert( std::is_same< typename _Vector1::value_type
	                           , typename _Vector2::value_type >::value
				 , "Element types of vectors must match!" );
	assert( is_row(X) or is_column(X) );
	assert( is_row(Y) or is_column(Y) );
	assert(    is_row(X) ? X.width() : X.height()
	        == is_row(Y) ? Y.width() : Y.height() );

	return tcm::import::dot
		( is_row(X) ? X.width() : X.height()
		, X.data(), is_row(X) ? X.ldim() : 1
		, Y.data(), is_row(Y) ? Y.ldim() : 1 );
}




// ============================================================================
// ||                                                                        ||
// ||                        MATRIX-VECTOR MULTIPLICATION                    ||
// ||                                                                        ||
// ============================================================================


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculates matrix-vector product.

/// This function if an abstraction from the ?GEMV functions in BLAS level 2.
/// It performs compile-time dispatching and calls SGEMV, DGEMV, CGEMV or
/// ZGEMV depending on the template parameters.
/// The following operation is thus performed:
/// \f[ Y := \alpha \cdot \mathcal{O}(A) \cdot X + \beta\cdot Y.\f]
///
/// \tparam _F        Real/complex field.
/// \tparam _Matrix   Matrix of type \p _F, that is 
///                   `typename _Matrix::value_type` must be well formed and
///                   return \p _F.
/// \tparam _Vector1  Vector of type \p _F.
/// \tparam _Vector2  Vector of type \p _F.
///
/// The following is assumed about _Matrix, _Vector1 and _Vector2 types:
/// 1) It has a value_type typedef, that is
///    `typename _Vector?::value_type` represents a type.
/// 2) is_row() and is_column() are defined for this type.
/// 3) `_Vector?::width()` returns the width of the vector if seen as a
///    matrix.
/// 4) `_Vector?::height()` returns the height of the vector if seen as a
///    matrix.
/// 5) `_Vector?::ldim()` returns the leading dimension of the vector.
/// 6) `_Vector?::data()` returns a pointer to the underlying array of
///    elements.
///
/// \warning This function currently only works for column vectors. It
///          on my TODO list to fix this.
///
/// \param[in]  op_A  \f$ \mathcal{O} \f$, i.e. the operator applied to \p A.
/// \param[in]  alpha Coefficient \f$\alpha\f$.
/// \param[in]  A     Matrix \f$ A \f$.
/// \param[in]  X     Vector \f$ X \f$.
/// \param[in]  beta  Coefficient \f$\beta\f$.
/// \param[out] Y     Vector \f$ Y \f$.
/// \exception Nothrow guarantee.
///////////////////////////////////////////////////////////////////////////////
template<class _F, class _Matrix, class _Vector1, class _Vector2>
inline
auto gemv( Operator const op_A
         , _F const alpha, _Matrix const& A, _Vector1 const& X
	     , _F const beta,  _Vector2 & Y ) noexcept -> void
{
	MEASURE;
	static_assert( std::is_same<typename _Matrix::value_type, _F>::value
				 , "Element type of _Matrix must match _F!" );
	static_assert( std::is_same<typename _Vector1::value_type, _F>::value
				 , "Element type of _Vector1 must match _F!" );
	static_assert( std::is_same<typename _Vector2::value_type, _F>::value
				 , "Element type of _Vector2 must match _F!" );

	assert( is_column(X) and is_column(Y) );
	assert( op_A == Operator::None 
		? A.height() == Y.height() and A.width()  == X.height()
	    : A.width()  == Y.height() and A.height() == X.height() );

	gemv( op_A, A.height(), A.width()
	    , alpha, A.data(), A.ldim(), X.data(), 1
	    , beta, Y.data(), 1 );
}







// ============================================================================
// ||                                                                        ||
// ||                    MATRIX-MATRIX MULTIPLICATION                        ||
// ||                                                                        ||
// ============================================================================


///////////////////////////////////////////////////////////////////////////////
/// \brief Calculates matrix-matrix product.

/// This function is an abstraction from ?GEMM functions in BLAS level 3. It
/// performs compile-time dispatching on template arguments and calls one of
/// SGEMM, DGEMM, CGEMM, ZGEMM.
/// The following operation is thus performed:
/// \f[ C := \alpha\cdot\mathcal{O}_A(A)\cdot\mathcal{O}_B(B)+\beta\cdot C.\f]
///
/// \tparam _F        Real/complex field.
/// \tparam _Matrix1  Matrix of type _F.
/// \tparam _Matrix2  Matrix of type _F.
/// \tparam _Matrix3  Matrix of type _F.
///
/// The following is assumed about _Matrix1, _Matrix2 and _Matrix3 types:
/// 1) It has a value_type typedef, that is
///    `typename _Matrix?::value_type` represents a type.
/// 2) is_row() and is_column() are defined for this type.
/// 3) `_Matrix?::width()` returns the width of the matrix.
/// 4) `_Matrix?::height()` returns the height of the matrix.
/// 5) `_Matrix?::ldim()` returns the leading dimension of the matrix.
/// 6) `_Matrix?::data()` returns a pointer to the underlying array of
///    elements.
///
/// \param[in]  op_A  \f$\mathcal{O}_A\f$, i.e. the operator applied to 
///                   matrix \f$ A \f$.
/// \param[in]  op_B  \f$\mathcal{O}_B\f$, i.e. the operator applied to 
///                   matrix \f$ B \f$.
/// \param[in]  alpha Coefficient \f$ \alpha \f$.
/// \param[in]  A     Matrix \f$ A \f$.
/// \param[in]  B     Matrix \f$ B \f$.
/// \param[in]  beta  Coefficient \f$ \beta \f$.
/// \param[out] C     Matrix \f$ C \f$.
/// \exception Nothrow guarantee.
/////////////////////////////////////////////////////////////////////////////// 
template<class _F, class _Matrix1, class _Matrix2, class _Matrix3>
inline
auto gemm( Operator const op_A, Operator const op_B
         , _F const alpha, _Matrix1 const& A, _Matrix2 const& B
         , _F const beta, _Matrix3 & C ) noexcept -> void
{
	MEASURE;
	static_assert( std::is_same<typename _Matrix1::value_type, _F>::value
	             , "Element type of _Matrix1 must match _F." );
	static_assert( std::is_same<typename _Matrix2::value_type, _F>::value
	             , "Element type of _Matrix2 must match _F." );
	static_assert( std::is_same<typename _Matrix3::value_type, _F>::value
	             , "Element type of _Matrix3 must match _F." );

#	ifndef NDEBUG
	auto _A_height_after = A.height();
	auto _A_width_after  = A.width();
	if(op_A != Operator::None) std::swap(_A_height_after, _A_width_after);

	auto _B_height_after = B.height();
	auto _B_width_after  = B.width();
	if(op_B != Operator::None) std::swap(_B_height_after, _B_width_after);

	assert(C.height() == _A_height_after);
	assert(C.width()  == _B_width_after);
	assert(_A_width_after == _B_height_after);
#	endif	

	gemm( op_A, op_B
	    , C.height()
	    , C.width()
	    , op_A == Operator::None ? A.width() : A.height()
	    , alpha, A.data(), A.ldim(), B.data(), B.ldim()
	    , beta, C.data(), C.ldim()
	    );
}



} // namespace blas

} // namespace tcm

#endif // TCM_BLAS_HPP
