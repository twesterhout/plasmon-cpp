#ifndef TCM_ITERATOR_HPP
#define TCM_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>


namespace tcm {


template<class _Tp>
class const_blas_iterator 
	: public boost::iterator_facade
	    < const_blas_iterator<_Tp>
	    , _Tp const
	    , boost::random_access_traversal_tag
		, _Tp
	    > {

private:
	using _facade_type = 
		         boost::iterator_facade< const_blas_iterator<_Tp>
	                                   , _Tp const
	                                   , boost::random_access_traversal_tag
	                                   , _Tp
	                                   >;
public:
	
	using value_type        = typename _facade_type::value_type;
	using difference_type   = typename _facade_type::difference_type;
	using reference         = typename _facade_type::reference;
	using size_type         = std::size_t;

	constexpr
	const_blas_iterator() noexcept
		: _data{ nullptr }, _step{ 1 }
	{
	}

	constexpr
	const_blas_iterator( value_type const* const data
	                   , difference_type const step = 1) noexcept
	    : _data{ data }, _step{ step }
	{ 
		assert(step != 0); 
	}

private:
	friend boost::iterator_core_access;

	constexpr auto dereference() const noexcept -> reference 
	{ assert(_data);
	  return *_data; }
	
	constexpr auto increment() noexcept 
	{ assert(_data);
	  _data += _step; }

	constexpr auto decrement() noexcept 
	{ assert(_data);
	  _data -= _step; }

	constexpr auto advance(difference_type const n) noexcept 
	{ assert(_data);
	  _data += n * _step; }

	constexpr auto distance_to(const_blas_iterator const& other) 
		const noexcept -> difference_type
	{ assert(_step == other._step);
	  return (other._data - _data) / _step; }

	constexpr auto equal(const_blas_iterator const& other) 
		const noexcept -> bool
	{ assert(_step == other._step);
	  return _data == other._data; }

private:
	value_type const* _data;
	difference_type   _step;
};




template<class _Tp>
class blas_iterator 
	: public boost::iterator_facade
	    < blas_iterator<_Tp>
	    , _Tp
	    , boost::random_access_traversal_tag
	    > {

private:
	using _facade_type = 
		typename boost::iterator_facade< const_blas_iterator<_Tp>
	                                   , _Tp
	                                   , boost::random_access_traversal_tag
	                                   >;
public:

	using value_type        = typename _facade_type::value_type;
	using difference_type   = typename _facade_type::difference_type;
	using reference         = typename _facade_type::reference;
	using size_type         = std::size_t;

	constexpr 
	blas_iterator() noexcept
		: _data{ nullptr }, _step{ 1 }
	{
	}

	constexpr
	blas_iterator( value_type* const data
	             , difference_type const step = 1 ) noexcept
	    : _data{ data }, _step{ step }
	{ 
		assert(step != 0); 
	}

private:
	friend boost::iterator_core_access;

	constexpr auto dereference() const noexcept -> reference
	{ assert(_data);
	  return *_data; }

	constexpr auto increment() noexcept
	{ assert(_data);
	  _data += _step; }

	constexpr auto decrement() noexcept
	{ assert(_data);
	  _data -= _step; }

	constexpr auto advance(difference_type n) noexcept
	{ _data += n * _step; }

	constexpr auto distance_to(blas_iterator const& other) 
		const noexcept -> difference_type
	{ assert(_step == other._step);
	  return (other._data - _data) / _step; }

	constexpr auto equal(blas_iterator const& other) 
		const noexcept -> bool
	{ assert(_step == other._step); 
	  return _data == other._data; 
	}

	operator const_blas_iterator<value_type>() const noexcept
	{ return {_data, _step}; }

private:
	value_type       _data;
	difference_type  _step;
};




} // namespace tcm


#endif // TCM_ITERATOR_HPP
