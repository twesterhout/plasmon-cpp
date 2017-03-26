#ifndef TCM_LAPACK_WRAPPER_MKL_HPP
#define TCM_LAPACK_WRAPPER_MKL_HPP


#include <complex>

#include <utils.hpp>

#define REGISTER_LAPACK(general_f, s_f, d_f, c_f, z_f)                      \
	namespace {                                                             \
		template<class... Args>                                             \
                __attribute__((always_inline))                              \
		inline                                                              \
		auto general_f##_impl( utils::Type2Type<float>                      \
		                     , Args&&... args) noexcept                     \
		{ return s_f(std::forward<Args>(args)...); }                        \
		                                                                    \
		template<class... Args>                                             \
                __attribute__((always_inline))                              \
		inline                                                              \
		auto general_f##_impl( utils::Type2Type<double>                     \
		                     , Args&&... args) noexcept                     \
		{ return d_f(std::forward<Args>(args)...); }                        \
		                                                                    \
		template<class... Args>                                             \
                __attribute__((always_inline))                              \
		inline                                                              \
		auto general_f##_impl( utils::Type2Type< std::complex<float> >      \
		                     , Args&&... args) noexcept                     \
		{ return c_f(std::forward<Args>(args)...); }                        \
		                                                                    \
		template<class... Args>                                             \
                __attribute__((always_inline))                              \
		inline                                                              \
		auto general_f##_impl( utils::Type2Type< std::complex<double> >     \
		                     , Args&&... args) noexcept                     \
		{ return z_f(std::forward<Args>(args)...); }                        \
	}                                                                       \
                                                                            \
	template<class T, class... Args>                                        \
	__attribute__((always_inline))                                          \
	inline                                                                  \
	auto general_f(Args&&... args) noexcept                                 \
	{ return general_f##_impl( utils::Type2Type<T>{}                        \
	                         , std::forward<Args>(args)...);                \
	}   


#ifndef USING_INTEL_MKL
#	error "Need Intel MKL"
#endif




namespace tcm { 

namespace import {

#define MKL_Complex8  std::complex<float>
#define MKL_Complex16 std::complex<double>

#include <mkl.h>

REGISTER_LAPACK(heev, ssyev_, dsyev_, cheev_, zheev_)
REGISTER_LAPACK(heevr, ssyevr_, dsyevr_, cheevr_, zheevr_)
REGISTER_LAPACK(geev, sgeev_, dgeev_, cgeev_, zgeev_)


} // namespace import

} // namespace tcm


#undef REGISTER_LAPACK



#endif // TCM_LAPACK_WRAPPER_MKL_HPP
