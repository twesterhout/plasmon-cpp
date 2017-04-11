#ifndef TCM_LAPACK_WRAPPER_ATLAS_HPP
#define TCM_LAPACK_WRAPPER_ATLAS_HPP


#include <complex>

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
        __attribute__((always_inline))                                      \
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





#ifndef USING_ATLAS
#	error "Need ATLAS"
#endif


namespace tcm {

namespace import {


extern "C" {


//                   ===================
//                   |      ?SYEV      |
//                   ===================

void ssyev_
    ( const char* JOB, const char* UPLO, const int* N
    , float* A, const int* LDA, float* W
    , float* WORK, const int* LWORK
    , int* INFO );

void dsyev_
    ( const char* JOB, const char* UPLO, const int* N
    , double* A, const int* LDA, double* W
    , double* WORK, const int* LWORK
    , int* INFO );




//                   ===================
//                   |      ?SYEVR     |
//                   ===================

void ssyevr_
    ( char const* JOBZ, char const* RANGE, char const* UPLO, int const* N
    , float* A, int const* LDA
    , float const* VL, float const* VU, int const* IL, int const* IU
    , float const* ABSTOL, int* M
    , float* W, float* Z, int const* LDZ
    , int* ISUPPZ
    , float* WORK, int const* LWORK
    , int* IWORK, int const* LIWORK
    , int* INFO );

void dsyevr_
    ( char const* JOBZ, char const* RANGE, char const* UPLO, int const* N
    , double* A, int const* LDA
    , double const* VL, double const* VU, int const* IL, int const* IU
    , double const* ABSTOL, int* M
    , double* W, double* Z, int const* LDZ
    , int* ISUPPZ
    , double* WORK, int const* LWORK
    , int* IWORK, int const* LIWORK
    , int* INFO );




//                   ===================
//                   |      ?HEEV      |
//                   ===================

void cheev_
    ( const char* JOB, const char* UPLO, const int* N
    , std::complex<float>* A, const int* LDA, float* W
    , std::complex<float>* WORK, const int* LWORK, float* RWORK
    , int* INFO );

void zheev_
    ( const char* JOB, const char* UPLO, const int* N
    , std::complex<double>* A, const int* LDA, double* W
    , std::complex<double>* WORK, const int* LWORK, double* RWORK
    , int* INFO );




//                   ===================
//                   |      ?HEEVR     |
//                   ===================

void cheevr_
    ( char const* JOBZ, char const* RANGE, char const* UPLO, int const* N
    , std::complex<float>* A, int const* LDA
    , float const* VL, float const* VU, int const* IL, int const* IU
    , float const* ABSTOL, int* M
    , float* W, std::complex<float>* Z, int const* LDZ
    , int* ISUPPZ
    , std::complex<float>* WORK, int* LWORK
    , float* RWORK, int* LRWORK
    , int* IWORK, int* LIWORK
    , int* INFO );

void zheevr_
    ( char const* JOBZ, char const* RANGE, char const* UPLO, int const* N
    , std::complex<double>* A, int const* LDA
    , double const* VL, double const* VU, int const* IL, int const* IU
    , double const* ABSTOL, int* M
    , double* W, std::complex<double>* Z, int const* LDZ
    , int* ISUPPZ
    , std::complex<double>* WORK, int* LWORK
    , double* RWORK, int* LRWORK
    , int* IWORK, int* LIWORK
    , int* INFO );




//                   ===================
//                   |      ?GEEV      |
//                   ===================

void sgeev_
    ( char const* JOBVL, char const* JOBVR, int const* N
    , float* A, int const* LDA, float* WR, float* WI
    , float* VL, int const* LDVL, float* VR, int const* LDVR
    , float* WORK, int const* LWORK
    , int* info );

void dgeev_
    ( char const* JOBVL, char const* JOBVR, int const* N
    , double* A, int const* LDA, double* WR, double* WI
    , double* VL, int const* LDVL, double* VR, int const* LDVR
    , double* WORK, int const* LWORK
    , int* info );

void cgeev_
    ( char const* JOBVL, char const* JOBVR, int const* N
    , std::complex<float>* A, int const* LDA, std::complex<float>* W
    , std::complex<float>* VL, int const* LDVL
    , std::complex<float>* VR, int const* LDVR
    , std::complex<float>* WORK, int const* LWORK, float* RWORK
    , int* info );

void zgeev_
    ( char const* JOBVL, char const* JOBVR, int const* N
    , std::complex<double>* A, int const* LDA, std::complex<double>* W
    , std::complex<double>* VL, int const* LDVL
    , std::complex<double>* VR, int const* LDVR
    , std::complex<double>* WORK, int const* LWORK, double* RWORK
    , int* info );

} // extern "C"



REGISTER_LAPACK(heev, ssyev_, dsyev_, cheev_, zheev_)
REGISTER_LAPACK(heevr, ssyevr_, dsyevr_, cheevr_, zheevr_)
REGISTER_LAPACK(geev, sgeev_, dgeev_, cgeev_, zgeev_)



} // namespace import

} // namespace tcm




#undef REGISTER_LAPACK

#endif // TCM_LAPACK_WRAPPER_ATLAS_HPP
