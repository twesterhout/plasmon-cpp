#ifndef TCM_MATRIX_HPP
#define TCM_MATRIX_HPP


#include <cassert>
#include <complex>
#include <memory>
#include <type_traits>
#include <algorithm>

#include <iterator.hpp>


///////////////////////////////////////////////////////////////////////////////
/// \file matrix.hpp
/// \brief This file implements a Matrix class and some operations associated
///        with it.
///////////////////////////////////////////////////////////////////////////////


namespace tcm {


namespace {
///////////////////////////////////////////////////////////////////////////////
/// \brief Contiguous chunk of memory. 

/// This class is like a %std::vector except that it allows no resizing and is
/// thus smaller in size by `sizeof(void*)`.
///
/// \tparam _Tp    element type.
/// \tparam _Alloc allocator type.
///////////////////////////////////////////////////////////////////////////////
template < class _Tp
         , class _Alloc
         >
struct _Storage 
	: public std::allocator_traits<_Alloc>::template rebind_alloc<_Tp> {

private:
	using _Tp_alloc_type    = typename std::allocator_traits<_Alloc>::template 
		rebind_alloc<_Tp>;
	using _Alloc_traits     = std::allocator_traits<_Tp_alloc_type>;

public:
	using value_type        = _Tp;
	using reference         = std::add_lvalue_reference_t<_Tp>;
	using const_reference   = std::add_const_t<reference>;
	using pointer           = typename _Alloc_traits::pointer;
	using const_pointer     = typename _Alloc_traits::const_pointer;
	using size_type         = typename _Alloc_traits::size_type;
	using difference_type   = typename _Alloc_traits::difference_type;
	using allocator_type    = _Alloc;

private:
	pointer _start;
	pointer _finish;

	auto _allocate(size_type const n) -> pointer
	{
		return n != 0 
			? _Alloc_traits::allocate(*this, n)
			: pointer{};
	}

	auto _deallocate(pointer const p, size_type const n) -> void
	{
		if (p) _Alloc_traits::deallocate(*this, p, n);
	}

	auto _create_storage(size_type const n)  -> void
	{
		_start  = _allocate(n);
		_finish = _start + n;
	}

	auto _swap_data(_Storage & x) noexcept -> void
	{
		using std::swap;
		swap(_start, x._start);
		swap(_finish, x._finish);
	}

	auto _get_Tp_allocator() noexcept -> _Tp_alloc_type &
	{
		return *static_cast<_Tp_alloc_type*>(this);
	}

	auto _get_Tp_allocator() const noexcept -> _Tp_alloc_type const&
	{
		return *static_cast<_Tp_alloc_type const*>(this);
	}


public:
	_Storage()
		noexcept( std::is_nothrow_default_constructible<_Tp_alloc_type>::value
		          and std::is_nothrow_constructible<pointer, nullptr_t>::value )
		: _Tp_alloc_type{}
		, _start{ nullptr }
		, _finish{ nullptr }
	{
	}

	_Storage(_Storage const& x)
		: _Storage{ x._get_Tp_allocator() }
	{
		_create_storage(x._finish - x._start);
		std::uninitialized_copy( x._start, x._finish
		                       , _start );
	}

	_Storage(allocator_type const& a)
		: _Tp_alloc_type{ a }
		, _start{ nullptr }
		, _finish{ nullptr }
	{
	}

	_Storage(size_type const n)
		: _Tp_alloc_type{}
	{
		_create_storage(n);
	}

	_Storage(allocator_type const& a, size_type const n)
		: _Tp_alloc_type{ a }
	{
		_create_storage(n);
	}

	_Storage(_Tp_alloc_type && a)
		: _Tp_alloc_type{ std::move(a) }
		, _start{ nullptr }
		, _finish{ nullptr }
	{
	}

	_Storage(_Storage && x)
		: _Storage{ std::move(x._get_Tp_allocator()) }
	{
		_swap_data(x);
	}

	~_Storage()
	{
		_deallocate(_start, _finish - _start);
	}

	constexpr auto size() const noexcept -> size_type
	{ return _finish - _start; }

	constexpr auto data() const noexcept -> pointer
	{ return _start; }

	friend auto swap(_Storage & x, _Storage & y)
	{
		std::swap(x._start, y._start);
		std::swap(x._finish, y._finish);
	}
};
} // unnamed namespace




///////////////////////////////////////////////////////////////////////////////
/// \brief A simple wrapper around the (data, ldim, width) representation
///        of matrices used in LAPACK.

/// This is a _container_ class in the sense that it manages its own memory.
/// Matrix class is meant to be used with LAPACK/BLAS and thus focuses on
/// fundamental numeric types.
///////////////////////////////////////////////////////////////////////////////
template< class _Tp
        , class _Alloc = std::allocator<_Tp>
        >
class Matrix {

private:
	using _Storage_type = _Storage<_Tp, _Alloc>;

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
	_Storage_type   _storage;
	size_type       _height;
	size_type       _width;
	difference_type _ldim;



public:

	///////////////////////////////////////////////////////////////////////////
	/// \brief Default constructor. 
	///////////////////////////////////////////////////////////////////////////
	Matrix()
		noexcept(std::is_nothrow_default_constructible<_Storage_type>::value)
	    : _storage{}
		, _height{ 0 }
		, _width{ 0 }
		, _ldim{ 0 }
	{
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Move constructor.
	///////////////////////////////////////////////////////////////////////////
	Matrix(Matrix && x)
		noexcept(std::is_nothrow_move_constructible<_Storage_type>::value)
	    : _storage{ std::move(x._storage) }
		, _height{ 0 }
		, _width{ 0 }
		, _ldim{ 0 }
	{
		std::swap(_height, x._height);
		std::swap(_width, x._width);
		std::swap(_ldim, x._ldim);
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Copy constructor.
	///////////////////////////////////////////////////////////////////////////
	Matrix(Matrix const& x)
		noexcept(std::is_nothrow_copy_constructible<_Storage_type>::value)
	    : _storage{ x._storage }
		, _height{ x._height }
		, _width{ x._width }
		, _ldim{ x._ldim }
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
		noexcept(std::is_nothrow_constructible< _Storage_type
		                                      , size_type>::value)
		: _storage{ height * width }
		, _height{ height }
		, _width{ width }
		, _ldim{ static_cast<difference_type>(height) }
	{
		assert( height < static_cast<size_type>(
			std::numeric_limits<difference_type>::max() ) );
	}

	///////////////////////////////////////////////////////////////////////////
	/// \brief Constructs a Matrix of given dimensions.

	/// \param height Height of the matrix, i.e. number of rows.
	/// \param width  Width of the matrix, i.e. number of columns.
	/// \param ldim   Leading dimension of the matrix. Must be non-zero.
	///
	/// \sa Matrix(size_type const, size_type const).
	///////////////////////////////////////////////////////////////////////////
	Matrix( size_type const height
	      , size_type const width
	      , difference_type const ldim )
		: Matrix{ height, width }
	{
		if (ldim == 0) 
			throw std::invalid_argument{"Leading dimension mustn't be zero."};
		_ldim = ldim;
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
	constexpr auto ldim() const noexcept -> difference_type
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
	{ return data() + (i + _ldim * j); }


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a non-const pointer to the element at row \p i and 
	/// column \p j.

	/// \se data(std::size_t const, std::size_t const) const
	///////////////////////////////////////////////////////////////////////////
	constexpr
	auto data( size_type const i
	         , size_type const j ) noexcept -> pointer
	{ return data() + (i + _height * j); }


	///////////////////////////////////////////////////////////////////////////
	/// Returns the element at row \p i and column \p j.

	/// \warning Mind you: not a const-reference, by the element by value!
	/// This is due to the fact that this class is meant to be used with
	/// elementary number types such as `float` or `int`.
	///////////////////////////////////////////////////////////////////////////
	constexpr
	auto operator() ( size_type const i
	                , size_type const j ) const noexcept -> value_type
	{ return *data(i, j); 
	}


	///////////////////////////////////////////////////////////////////////////
	/// Returns a non-const reference to the element at row \p i and 
	/// column \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr
	auto operator() ( size_type const i
	                , size_type const j ) noexcept -> reference
	{ return *data(i, j);
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant row iterator to the beginning of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cbegin_row(size_type const i) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(i < _height);
	  return const_blas_iterator<value_type>{data(i, 0), ldim()};
	}
	
	
	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating row iterator to the beginning of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto begin_row(size_type const i) noexcept
		-> blas_iterator<value_type>
	{ assert(i < _height);
	  return blas_iterator<value_type>{data(i, 0), ldim()};
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant row iterator to the end of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cend_row(size_type const i) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(i < _height);
	  return const_blas_iterator<value_type>{data(i, width()), ldim()};
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating row iterator to the end of row \p i.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto end_row(size_type const i) noexcept
		-> blas_iterator<value_type>
	{ assert(i < _height);
	  return blas_iterator<value_type>{data(i, width()), ldim()};
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant column iterator to the beginning of column 
	/// \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cbegin_column(size_type const j) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(j < _width);
	  return {data(0, j), 1};
	}
	

	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating column iterator to the beginning of column 
	/// \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto begin_column(size_type const j) noexcept
		-> blas_iterator<value_type>
	{ assert(j < _width);
	  return blas_iterator<value_type>{data(0, j), 1};
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a constant column iterator to the end of column \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto cend_column(size_type const j) const noexcept
		-> const_blas_iterator<value_type>
	{ assert(j < _width);
	  return const_blas_iterator<value_type>{data(height(), j), 1};
	}


	///////////////////////////////////////////////////////////////////////////
	/// \brief Returns a mutating column iterator to the end of column \p j.
	///////////////////////////////////////////////////////////////////////////
	constexpr auto end_column(size_type const j) noexcept
		-> blas_iterator<value_type>
	{ assert(j < _width);
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


template <class _T, class _Alloc>
auto operator<< (std::ostream& out, Matrix<_T, _Alloc> const& A) 
	-> std::ostream&
{
	for (std::size_t i = 0; i < A.height(); ++i) {
		for (std::size_t j = 0; j < A.width(); ++j)
			out << A(i, j) << '\t';
		out << '\n';
	}
	return out;
}


template<class _T, class _Alloc>
auto operator>> (std::istream& in, Matrix<_T, _Alloc>& A) 
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
