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
	using boost::serialization::make_array;

	if (version == 0) {
		auto const height = matrix.height();
		auto const width  = matrix.width();
		ar << height;
		ar << width;
		ar << make_array(matrix.data(), matrix.height() * matrix.width());
	}
	else {
		auto const height = matrix.height();
		auto const width  = matrix.width();
		auto const ldim   = matrix.ldim();
		ar << height;
		ar << width;
		ar << ldim;
		
		for (std::size_t col = 0; col < matrix.width(); ++col) {
			for (std::size_t row = 0; row < matrix.height(); ++row) {
				// we need a reference here
				ar << *matrix.data(row, col);
			}
		}
	}
}


template <class _Archive, class _Tp, class _Alloc>
auto load( _Archive & ar
         , tcm::Matrix<_Tp, _Alloc> & matrix
         , unsigned int version )
{
	using M = tcm::Matrix<_Tp, _Alloc>;
	using boost::serialization::make_array;

	if (version == 0) {
		typename M::size_type height, width;

		ar >> height;
		ar >> width;
		M temp{height, width};
		ar >> make_array(temp.data(), height * width);
		swap(matrix, temp);
	}
	else {
		typename M::size_type       height, width;
		typename M::difference_type ldim;
		ar >> height;
		ar >> width;
		ar >> ldim;
	
		M temp{height, width, ldim};
		for (std::size_t col = 0; col < matrix.width(); ++col) {
			for (std::size_t row = 0; row < matrix.height(); ++row) {
				ar >> temp(row, col);
			}
		}
		swap(matrix, temp);
	}
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
