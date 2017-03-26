#ifndef TCM_HERMITIAN_HPP
#define TCM_HERMITIAN_HPP

#include <cassert>
#include <limits>
#include <memory>


#include <lapack_wrapper.hpp>
#include <benchmark.hpp>




// ============================================================================
// ||                                                                        ||
// ||                             ? H E E V R                                ||
// ||                                                                        ||
// ============================================================================


namespace tcm {

namespace lapack {


namespace {


template< class _T
        , class = std::enable_if_t
                  <    std::is_same<_T, float>()
                    or std::is_same<_T, double>()
		          >
        >
auto heevr_impl( int const N
               , _T* A, int const LDA
               , _T* W
               , _T* Z, int const LDZ
               , utils::Type2Type<_T> ) -> void
{
	if (N == 0) return;

	assert(N > 0);
	assert(A != nullptr and LDA >= N);
	assert(W != nullptr);
	assert(LDZ >= (Z == nullptr) ? 1 : N);

	char const JOBZ   = (Z == nullptr) ? 'N' : 'V';
	char const RANGE  = 'A';
	char const UPLO   = 'U';
	_T   const VL     { 0.0 };
	_T   const VU     { 0.0 };
	int  const IL     = 0;
	int  const IU     = 0;
	_T   const ABSTOL { 0.0 }; //std::numeric_limits<_T>::min();
	int        M      = 0;
	auto       ISUPPZ = std::make_unique<int[]>(2 * N);
	int        LWORK  = -1;
	int        LIWORK = -1;
	int        INFO   = 0;

	{
		_T  _work_dummy;
		int _iwork_dummy;

		tcm::import::heevr<_T>
			( &JOBZ, &RANGE, &UPLO, &N
			, A, &LDA
			, &VL, &VU
			, &IL, &IU
			, &ABSTOL, &M
			, W
			, Z, &LDZ
			, ISUPPZ.get()
			, &_work_dummy, &LWORK
			, &_iwork_dummy, &LIWORK
			, &INFO
		    );

		LWORK = static_cast<int>(_work_dummy);
		LIWORK = _iwork_dummy;
	}

	auto       WORK   = std::make_unique<_T[]>(LWORK);
	auto       IWORK  = std::make_unique<int[]>(LIWORK);

	tcm::import::heevr<_T>
		( &JOBZ, &RANGE, &UPLO, &N
		, A, &LDA
		, &VL, &VU
		, &IL, &IU
		, &ABSTOL, &M
		, W, Z, &LDZ
		, ISUPPZ.get()
		, WORK.get(), &LWORK
		, IWORK.get(), &LIWORK
		, &INFO
		);

	if (INFO < 0) {
		throw std::invalid_argument{ "Argument #" + std::to_string(-INFO) 
		                           + " had an illegal value." };
	} 
	else if (INFO > 0) {
		throw std::runtime_error{"Call to ?SYEVR failed."};
	}
}



template< class _T
        , class = std::enable_if_t
                  <    std::is_same<_T, float>()
                    or std::is_same<_T, double>()
                  >
        >
auto heevr_impl( int const N
               , std::complex<_T>* A, int const LDA
               , _T* W
               , std::complex<_T>* Z, int const LDZ
               , utils::Type2Type<std::complex<_T>> ) -> void
{
	if (N == 0) return;

	assert(N > 0);
	assert(A != nullptr and LDA >= N);
	assert(W != nullptr);
	assert(LDZ >= (Z == nullptr) ? 1 : N);

	char const JOBZ   = (Z == nullptr) ? 'N' : 'V';
	char const RANGE  = 'A';
	char const UPLO   = 'U';
	_T   const VL     { 0.0 };
	_T   const VU     { 0.0 };
	int  const IL     = 0;
	int  const IU     = 0;
	_T   const ABSTOL { 0.0 }; //std::numeric_limits<T>::min();
	int        M      = 0;
	auto       ISUPPZ = std::make_unique<int[]>(2 * N);
	int        LWORK  = -1;
	int        LRWORK = -1;
	int        LIWORK = -1;
	int        INFO   = 0;

	{
		std::complex<_T>   _work_dummy;
		_T                 _rwork_dummy;
		int                _iwork_dummy;

		tcm::import::heevr<std::complex<_T>>
			( &JOBZ, &RANGE, &UPLO, &N
			, A, &LDA
			, &VL, &VU
			, &IL, &IU
			, &ABSTOL, &M
			, W
			, Z, &LDZ
			, ISUPPZ.get()
			, &_work_dummy, &LWORK
			, &_rwork_dummy, &LRWORK
			, &_iwork_dummy, &LIWORK
			, &INFO
		    );

		LWORK  = static_cast<int>(std::real(_work_dummy));
		LRWORK = static_cast<int>(_rwork_dummy);
		LIWORK = _iwork_dummy;
	}

	auto       WORK   = std::make_unique<std::complex<_T>[]>(LWORK);
	auto       RWORK  = std::make_unique<_T[]>(LRWORK);
	auto       IWORK  = std::make_unique<int[]>(LIWORK);

	tcm::import::heevr<std::complex<_T>>
		( &JOBZ, &RANGE, &UPLO, &N
		, A, &LDA
		, &VL, &VU
		, &IL, &IU
		, &ABSTOL, &M
		, W, Z, &LDZ
		, ISUPPZ.get()
		, WORK.get(), &LWORK
		, RWORK.get(), &LRWORK
		, IWORK.get(), &LIWORK
		, &INFO
		);

	if(INFO < 0) {
		throw std::invalid_argument{ "Argument #" + std::to_string(-INFO) 
		                           + " had an illegal value." };
	} 
	else if(INFO > 0) {
		throw std::runtime_error{"Call to ?HEEVR failed."};
	}
}

} // unnamed namespace 



template<class _T>
auto heevr( std::size_t const n
          , _T* A, std::ptrdiff_t const lda
          , utils::Base<_T>* W
          , _T* Z, std::ptrdiff_t const ldz ) -> void
{
	MEASURE;

	heevr_impl( boost::numeric_cast<int>(n)
	          , A, boost::numeric_cast<int>(lda)
	          , W
	          , Z, boost::numeric_cast<int>(ldz)
	          , utils::Type2Type<_T>{} );
}



} // namespace lapack

} // namespace tcm








// ============================================================================
// ||                                                                        ||
// ||                              ? H E E V                                 ||
// ||                                                                        ||
// ============================================================================



namespace tcm {

namespace lapack {


namespace {


template< class _T
        , class = std::enable_if_t
                  <   std::is_same<_T, float>() 
                   or std::is_same<_T, double>()
                  >
        >
auto heev_impl( int const N
              , _T* A, int const LDA
              , _T* W
              , bool compute_eigenvectors
              , utils::Type2Type<_T> ) -> void
{
	if (N == 0) return;

	assert(N > 0);
	assert(A != nullptr and LDA >= N);
	assert(W != nullptr);

	char const JOBZ  = compute_eigenvectors ? 'V' : 'N';
	char const UPLO  = 'U';
	int        LWORK = -1;
	int        INFO;

	{
		_T _work_dummy;

		tcm::import::heev<_T>
			( &JOBZ, &UPLO, &N
			, A, &LDA
			, W
			, &_work_dummy,  &LWORK
			, &INFO
			);

		LWORK = static_cast<int>(_work_dummy);
	}

	auto       WORK  = std::make_unique<_T[]>(LWORK);

	tcm::import::heev<_T>
		( &JOBZ, &UPLO, &N
		, A, &LDA
		, W
		, WORK.get(), &LWORK
		, &INFO
		);

	if (INFO < 0) {
		throw std::invalid_argument{ "Argument #" + std::to_string(-INFO) 
		                           + " had an illegal value." };
	}
	else if (INFO > 0) {
		throw std::runtime_error{"Call to ?SYEV failed."};
	}
}


template< class _T
        , class = std::enable_if_t
                  <   std::is_same<_T, float>() 
                   or std::is_same<_T, double>()
                  >
        >
auto heev_impl( int const N
              , std::complex<_T>* A, int const LDA
              , _T* W
              , bool compute_eigenvectors
              , utils::Type2Type<std::complex<_T>> ) -> void
{
	if(N == 0) return;

	assert(N > 0);
	assert(A != nullptr and LDA >= N);
	assert(W != nullptr);

	char const JOBZ  = compute_eigenvectors ? 'V' : 'N';
	char const UPLO  = 'U';
	auto       RWORK = std::make_unique<_T[]>(3 * N - 2);
	int        LWORK = -1;
	int        INFO;

	{
		std::complex<_T> _work_dummy;

		tcm::import::heev<std::complex<_T>>
			( &JOBZ, &UPLO, &N
			, A, &LDA
			, W
			, &_work_dummy, &LWORK
			, RWORK.get()
			, &INFO
			);

		LWORK = static_cast<int>(std::real(_work_dummy));
	}
	
	auto       WORK  = std::make_unique<std::complex<_T>[]>(LWORK);

	tcm::import::heev<std::complex<_T>>
		( &JOBZ, &UPLO, &N
		, A, &LDA
		, W
		, WORK.get(),  &LWORK
		, RWORK.get()
		, &INFO
		);

	if (INFO < 0) {
		throw std::invalid_argument{ "Argument #" + std::to_string(-INFO) 
		                           + " had an illegal value." };
	}
	else if (INFO > 0) {
		throw std::runtime_error{"Call to ?HEEV failed."};
	}
}

} // unnamed namespace


template<class _T>
inline
auto heev( std::size_t const n
         , _T* A, std::ptrdiff_t const lda
         , utils::Base<_T>* W
         , bool compute_eigenvectors ) -> void
{
	MEASURE;

	heev_impl( boost::numeric_cast<int>(n)
	         , A, boost::numeric_cast<int>(lda)
	         , W
	         , compute_eigenvectors
	         , utils::Type2Type<_T>{} );
}



} // namespace lapack

} // namespace tcm




#endif // TCM_HERMITIAN_HPP
