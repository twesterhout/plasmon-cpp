#ifndef TCM_MATRIX_HPP
#define TCM_MATRIX_HPP


#include <cassert>
#include <complex>
#include <memory>
#include <type_traits>
#include <algorithm>

#include <boost/align/aligned_allocator_adaptor.hpp>

#include <detail/utils.hpp>
#include <detail/iterator.hpp>


///////////////////////////////////////////////////////////////////////////////
/// \file include/matrix.hpp
/// \brief This file implements a Matrix class and some operations associated
///        with it.
///////////////////////////////////////////////////////////////////////////////


namespace tcm {


///////////////////////////////////////////////////////////////////////////////
/// \brief A simple wrapper around the (data, ldim, width) representation
///        of matrices used in LAPACK. Uses column-major ordering.

/// This is a _container_ class in the sense that it manages its own memory.
/// Matrix class is meant to be used with LAPACK/BLAS and thus focuses on
/// fundamental numeric types. There are two important differences
/// between this class and, for example, an %std::vector.
/// 1) As we plan to use matrices with LAPACK/BLAS routines, it helps
///    (from the performance point of view) to correctly align our matrices.
///    From Intel MKL's manual:
///    > To improve performance of your application that calls Intel MKL, 
///    > align your arrays on 64-byte boundaries and ensure that the leading 
///    > dimensions of the arrays are divisible by 64/element_size, where 
///    > element_size is the number of bytes for the matrix elements 
///    > (4 for single-precision real, 8 for double-precision real and 
///    > single-precision complex, and 16 for double-precision complex). 
/// 2) We don't initialise storage.
///
/// \tparam _Tp Element type of the matrix.
/// \tparam _Align Alignment. Must be a power of 2. This argument is used 
///                for two things:
///                1) Allocation of storage that is aligned at least to _Alloc.
///                2) Computing the leading dimension of the matrix to ensure
///                   that `this->data(0, 1)` is again aligned to _Alloc.
/// 
///                We use %boost::alignment::aligned_allocator_adaptor
///                for this.
///
///////////////////////////////////////////////////////////////////////////////
template< class _Tp
        , std::size_t _Align = 64 // std::alignment_of<_Tp>::value
        , class _Alloc = std::allocator<_Tp>
        >
class Matrix {

private:
	using _Aligned_Alloc = boost::alignment::
		aligned_allocator_adaptor<_Alloc, _Align>;
	using _Storage_type = utils::_Storage<_Tp, _Aligned_Alloc>;

	static_assert( (_Align != 0) and ((_Align & (_Align - 1)) == 0)
	             , "_Align must be a power of 2" );
	static_assert( _Align % std::alignment_of<_Tp>::value == 0
	             , "Invalid alignment." );

public:
	using value_type        = _Tp;
	using reference         = typename _Storage_type::reference;
	using const_reference   = typename _Storage_type::const_reference;
	using pointer           = typename _Storage_type::pointer;
	using const_pointer     = typename _Storage_type::const_pointer;
	using size_type         = typename _Storage_type::size_type;
	using difference_type   = typename _Storage_type::difference_type;
	using allocator_type    = typename _Storage_type::allocator_type;

private:
	size_type       _height;
	size_type       _width;
	size_type       _ldim;
	_Storage_type   _storage;

	static constexpr auto round_up(size_type const n) -> size_type
	{
		constexpr auto multiple = _Align / std::alignment_of<_Tp>::value;
		return ((n + multiple - 1) / multiple) * multiple;
	}

public:

	///////////////////////////////////////////////////////////////////////////
	/// \brief Default constructor. 
	///////////////////////////////////////////////////////////////////////////
	Matrix()
		noexcept(std::is_nothrow_default_constructible<_Storage_type>::value)
		: _height{ 0 }
		, _width{ 0 }
		, _ldim{ 0 }
	    , _storage{}
	{
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Move constructor.
	///////////////////////////////////////////////////////////////////////////
	Matrix(Matrix && x)
		noexcept(std::is_nothrow_move_constructible<_Storage_type>::value)
		: _height{ 0 }
		, _width{ 0 }
		, _ldim{ 0 }
	    , _storage{ std::move(x._storage) }
	{
		using std::swap;
		swap(_height, x._height);
		swap(_width, x._width);
		swap(_ldim, x._ldim);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Copy constructor.
	///////////////////////////////////////////////////////////////////////////
	Matrix(Matrix const& x)
		noexcept(std::is_nothrow_copy_constructible<_Storage_type>::value)
		: _height{ x._height }
		, _width{ x._width }
		, _ldim{ x._ldim }
	    , _storage{ x._storage }
	{
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Constructs a Matrix of given dimensions.

	/// \param height Height of the matrix, i.e. number of rows.
	/// \param width  Width of the matrix, i.e. number of columns.
	///
	/// \sa Matrix(size_type const, size_type const, difference_type const).
	///////////////////////////////////////////////////////////////////////////
	Matrix(size_type const height, size_type const width)
		: _height{ height }
		, _width{ width }
		// , _ldim{ height }
		, _ldim{ round_up(height) }
		, _storage{ round_up(height) * width }
	{
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Destructor.
	///////////////////////////////////////////////////////////////////////////
	~Matrix()
	{
	}


	///////////////////////////////////////////////////////////////////////////
	/// Move assignment operator. 
	
	/// No copy assignment operator is provided
	/// to make it less convenient to perform copies all the time.
	///////////////////////////////////////////////////////////////////////////
	auto operator=(Matrix && other) noexcept -> Matrix &
	{ swap(*this, other);
	  return *this;
	}

	///////////////////////////////////////////////////////////////////////////
	/// Returns the height of the matrix.

	/// \sa width(), ldim().
	///////////////////////////////////////////////////////////////////////////
	constexpr auto height() const noexcept -> size_type
	{ return _height; }


	///////////////////////////////////////////////////////////////////////////
	/// Returns the width of the matrix.
	
	/// \sa height(), ldim().
	///////////////////////////////////////////////////////////////////////////
	constexpr auto width() const noexcept -> size_type
	{ return _width; }


	///////////////////////////////////////////////////////////////////////////
	/// Returns the leading dimension of the matrix.

	/// \sa width(), height().
	///////////////////////////////////////////////////////////////////////////
	constexpr auto ldim() const noexcept -> size_type
	{ return _ldim; }


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a const pointer to the underlying one-dimensional array
	/// of elements.

	/// \sa data()
	///////////////////////////////////////////////////////////////////////////
	constexpr auto data() const noexcept -> const_pointer
	{ return _storage.data(); }


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a non-const pointer to the underlying one-dimensional 
	/// array of elements.

	/// \sa data() const
	///////////////////////////////////////////////////////////////////////////
	constexpr auto data() noexcept -> pointer
	{ return _storage.data(); }


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a const pointer to the element at row \p i and 
	/// column \p j.

	/// \se data(std::size_t const, std::size_t const)
	///////////////////////////////////////////////////////////////////////////
	constexpr
	auto data( size_type const i
	         , size_type const j ) const noexcept -> const_pointer
	{ return data() + (i + ldim() * j); }


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a non-const pointer to the element at row \p i and 
	/// column \p j.

	/// \se data(std::size_t const, std::size_t const) const
	///////////////////////////////////////////////////////////////////////////
	constexpr
	auto data( size_type const i
	         , size_type const j ) noexcept -> pointer
	{ return data() + (i + ldim() * j); }


	///////////////////////////////////////////////////////////////////////////
	/// Returns the element at row \p i and column \p j.

	/// \warning Mind you: not a const-reference, by the element by value!
	/// This is due to the fact that this class is meant to be used with
	/// elementary number types such as `float` or `int`.
	///////////////////////////////////////////////////////////////////////////
	constexpr
	auto operator() ( size_type const i
	                , size_type const j ) const noexcept -> value_type
	{ assert(i < height() and j < width() and "Index out of bounds.");
	  return *data(i, j); 
	}


	///////////////////////////////////////////////////////////////////////////
	/// Returns a non-const reference to the element at row \p i and 
	/// column \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr
	auto operator() ( size_type const i
	                , size_type const j ) noexcept -> reference
	{ assert(i < height() and j < width() and "Index out of bounds.");
	  return *data(i, j);
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant row iterator to the beginning of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cbegin_row(size_type const i) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(i < height() and "Index out of bounds.");
	  return const_blas_iterator<value_type>{ data(i, 0),
	      static_cast<difference_type>(ldim()) };
	}
	
	
	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating row iterator to the beginning of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto begin_row(size_type const i) noexcept
		-> blas_iterator<value_type>
	{ assert(i < height() and "Index out of bounds.");
	  return blas_iterator<value_type>{ data(i, 0), 
	      static_cast<difference_type>(ldim()) };
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant row iterator to the end of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cend_row(size_type const i) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(i < height() and "Index out of bounds.");
	  return const_blas_iterator<value_type>{ data(i, width()),
	      static_cast<difference_type>(ldim()) };
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating row iterator to the end of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto end_row(size_type const i) noexcept
		-> blas_iterator<value_type>
	{ assert(i < height() and "Index out of bounds.");
	  return blas_iterator<value_type>{ data(i, width()), 
	      static_cast<difference_type>(ldim()) };
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant column iterator to the beginning of column 
	/// \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cbegin_column(size_type const j) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(j < width() and "Index out of bounds.");
	  return {data(0, j), 1};
	}
	

	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating column iterator to the beginning of column 
	/// \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto begin_column(size_type const j) noexcept
		-> blas_iterator<value_type>
	{ assert(j < width() and "Index out of bounds.");
	  return blas_iterator<value_type>{data(0, j), 1};
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant column iterator to the end of column \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cend_column(size_type const j) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(j < width() and "Index out of bounds.");
	  return const_blas_iterator<value_type>{data(height(), j), 1};
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating column iterator to the end of column \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto end_column(size_type const j) noexcept
		-> blas_iterator<value_type>
	{ assert(j < width() and "Index out of bounds.");
	  return blas_iterator<value_type>{data(height(), j), 1};
	}

	
private:
	friend
	auto swap(Matrix<value_type>& lhs, Matrix<value_type>& rhs) noexcept
	{ using std::swap;
	  swap(lhs._storage, rhs._storage);
	  swap(lhs._height, rhs._height);
	  swap(lhs._width, rhs._width);
	  swap(lhs._ldim, rhs._ldim);
	}
};



///////////////////////////////////////////////////////////////////////////////
/// Returns true if the given matrix is just a row vector.
///////////////////////////////////////////////////////////////////////////////
template <class _M>
constexpr 
auto is_row(_M const& A) noexcept -> bool { return A.height() == 1; }

///////////////////////////////////////////////////////////////////////////////
/// Returns true if the given matrix is just a column vector.
///////////////////////////////////////////////////////////////////////////////
template <class _M>
constexpr 
auto is_column(_M const& A) noexcept -> bool { return A.width() == 1; }


///////////////////////////////////////////////////////////////////////////////
/// Returns true if the given matrix is a square matrix.
///////////////////////////////////////////////////////////////////////////////
template <class _M>
constexpr 
auto is_square(_M const& A) noexcept -> bool {return A.width() == A.height();}


template <class _T, std::size_t _Align, class _Alloc>
auto operator<< (std::ostream& out, Matrix<_T, _Align, _Alloc> const& A) 
	-> std::ostream&
{
	for (std::size_t i = 0; i < A.height(); ++i) {
		for (std::size_t j = 0; j < A.width(); ++j)
			out << A(i, j) << '\t';
		out << '\n';
	}
	return out;
}


template<class _T, std::size_t _Align, class _Alloc>
auto operator>> (std::istream& in, Matrix<_T, _Align, _Alloc>& A) 
	-> std::istream&
{
	for (std::size_t i = 0; i < A.height(); ++i) {
		for (std::size_t j = 0; j < A.width(); ++j) {
			in >> A(i, j);
		}
	}
	return in;
}


template<class Func>
auto build_matrix( std::size_t const height, std::size_t const width
                 , Func f )
{
	using T = decltype( f( std::declval<std::size_t>()
	                     , std::declval<std::size_t>()
	                     ) );
	Matrix<T> result{height, width};
	for (std::size_t j = 0; j < width; ++j)
		for (std::size_t i = 0; i < height; ++i)
			result(i, j) = f(i, j);
	return result;
}

} // namespace tcm




#endif // TCM_MATRIX_HPP
