#ifndef TCM_UTILS_HPP
#define TCM_UTILS_HPP


#include <iostream>
#include <type_traits>
#include <memory>
#include <utility>

namespace tcm {


namespace utils {


///////////////////////////////////////////////////////////////////////////////
/// \brief Contiguous chunk of memory. 

/// This class is like a %std::vector except that it allows no resizing and is
/// thus smaller in size by `sizeof(pointer)`.
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

private:
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
		// Case n == 0 is handled inside _allocate(),
		// no need to do it again here.
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
		std::copy( x._start, x._finish, _start );
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





template<class Enum>
constexpr
auto to_integral(Enum x) noexcept
{ return static_cast<std::underlying_type_t<Enum>>(x); }


template<class... Types> struct Type2Type {};


namespace {
template<class T> struct Base_impl { using type = T; };
template<class T> struct Base_impl<std::complex<T>> { using type = T; };
} // unnamed namespace


template<class T>
using Base = typename Base_impl<T>::type;



template <class _Alloc>
auto allocate_workspace(_Alloc& alloc, std::size_t n) 
{
	assert(n != 0);

	using _Alloc_traits = std::allocator_traits<_Alloc>;
	auto deleter = [&alloc, n] (auto* const p) {
		_Alloc_traits::deallocate(alloc, p, n);
	};

	return std::unique_ptr< typename _Alloc_traits::value_type
	                      , decltype(deleter) >
		( _Alloc_traits::allocate(alloc, n)
		, deleter );
}






template<class _F>
auto assert_valid(_F const x) noexcept -> void
{
#	ifndef NO_CHECK_NANS
	assert(not std::isnan(std::real(x)) and not std::isnan(std::imag(x)));
	assert(not std::isinf(std::real(x)) and not std::isinf(std::imag(x)));
#	endif
}


template<class _Iterator>
auto assert_valid(_Iterator begin, _Iterator end) -> void
{
#	ifndef NO_CHECK_NANS
	assert( not std::count_if
		( begin, end
	    , [] (auto const y) { 
		      return std::isnan(std::real(y))    
	              or std::isnan(std::imag(y))
				  or std::isinf(std::real(y))
				  or std::isinf(std::imag(y)); 
		  }
		) 
	);
#	endif
}


} // namespace utils







namespace hacking {


namespace {

template<class F>
struct operator_wrapper {
	F f;
};

template<class Lhs, class F>
struct operator_proxy {
	operator_wrapper<F> op;
	Lhs lhs;
};

template<class Lhs, class F>
auto operator< (Lhs&& lhs, operator_wrapper<F> op)
{ return operator_proxy<Lhs, F>{op, std::forward<Lhs>(lhs)}; }

template<class Lhs, class F, class Rhs>
auto operator> (operator_proxy<Lhs, F>&& proxy, Rhs&& rhs)
{ return proxy.op.f(proxy.lhs, std::forward<Rhs>(rhs)); }

} // unnamed namespace

template<class F>
auto register_operator(F f)
{ return operator_wrapper<F>{f}; }


} // namespace hacking








} // namespace tcm

#endif // TCM_UTILS_HPP
