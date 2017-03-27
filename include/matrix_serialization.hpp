#ifndef TCM_MATRIC_SERIALIZATION_HPP
#define TCM_MATRIC_SERIALIZATION_HPP


#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/split_free.hpp>


#include <matrix.hpp>


namespace boost {
namespace serialization {

template <class _Archive, class _Tp, class _Alloc>
auto save( _Archive & ar
         , tcm::Matrix<_Tp, _Alloc> const& matrix
         , unsigned int version )
{
	auto const height = matrix.height();
	auto const width  = matrix.width();
	ar << height;
	ar << width;

	for (std::size_t col = 0; col < matrix.width(); ++col) {
		for (std::size_t row = 0; row < matrix.height(); ++row) {
			// we need a reference here
			ar << *matrix.data(row, col);
		}
	}
}


template <class _Archive, class _Tp, class _Alloc>
auto load( _Archive & ar
         , tcm::Matrix<_Tp, _Alloc> & matrix
         , unsigned int version )
{
	using std::swap;
	using M = tcm::Matrix<_Tp, _Alloc>;

	typename M::size_type height, width;
	ar >> height;
	ar >> width;

	M temp{height, width};
	for (std::size_t col = 0; col < width; ++col) {
		for (std::size_t row = 0; row < height; ++row) {
			ar >> temp(row, col);
		}
	}
	swap(matrix, temp);
}


template<class _Archive, class _Tp, class _Alloc>
inline 
auto serialize( _Archive & ar
              , tcm::Matrix<_Tp, _Alloc> & matrix
              , unsigned int const file_version ) -> void
{
	split_free(ar, matrix, file_version); 
}


} // namespace serialization
} // namespace boost




#endif // TCM_MATRIC_SERIALIZATION_HPP
