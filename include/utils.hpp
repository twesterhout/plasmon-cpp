#ifndef TCM_UTILS_HPP
#define TCM_UTILS_HPP


#include <iostream>
#include <type_traits>
#include <memory>
#include <utility>

namespace tcm {


namespace utils {


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
