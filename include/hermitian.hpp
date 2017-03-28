#ifndef TCM_HERMITIAN_HPP
#define TCM_HERMITIAN_HPP

#include <cassert>
#include <limits>
#include <memory>

#include <boost/numeric/conversion/cast.hpp>

#include <utils.hpp>
#include <benchmark.hpp>
#include <lapack_wrapper.hpp>



// ============================================================================
// ||                                                                        ||
// ||                             ? H E E V R                                ||
// ||                                                                        ||
// ============================================================================


namespace tcm {

namespace lapack {


namespace {


template< class _Alloc
        , class _T
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

	/*
	using _Alloc_traits = std::allocator_traits<_Alloc>;
	typename _Alloc_traits::template rebind_alloc<_T>
		_T_alloc;
	typename _Alloc_traits::template rebind_alloc<int>
		_int_alloc;
	*/
	using size_type = std::make_unsigned_t<int>;

	char const JOBZ   = (Z == nullptr) ? 'N' : 'V';
	char const RANGE  = 'A';
	char const UPLO   = 'U';
	_T   const VL     { 0.0 };
	_T   const VU     { 0.0 };
	int  const IL     = 0;
	int  const IU     = 0;
	_T   const ABSTOL { 0.0 }; //std::numeric_limits<_T>::min();
	int        M      = 0;
	auto       ISUPPZ = 
		tcm::utils::_Storage<int, _Alloc>{2 * static_cast<size_type>(N)};
	// tcm::utils::allocate_workspace(_int_alloc, 2 * N);
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
			, ISUPPZ.data()
			, &_work_dummy, &LWORK
			, &_iwork_dummy, &LIWORK
			, &INFO
		    );

		LWORK = static_cast<int>(_work_dummy);
		LIWORK = _iwork_dummy;
	}

	auto       WORK   = 
		tcm::utils::_Storage<_T, _Alloc>{static_cast<size_type>(LWORK)};
	// tcm::utils::allocate_workspace(_T_alloc, LWORK);
	auto       IWORK  = 
		tcm::utils::_Storage<int, _Alloc>{static_cast<size_type>(LIWORK)};
	// allocate_workspace(_int_alloc, LIWORK);

	tcm::import::heevr<_T>
		( &JOBZ, &RANGE, &UPLO, &N
		, A, &LDA
		, &VL, &VU
		, &IL, &IU
		, &ABSTOL, &M
		, W, Z, &LDZ
		, ISUPPZ.data()
		, WORK.data(), &LWORK
		, IWORK.data(), &LIWORK
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



template< class _Alloc
        , class _T
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

	/*
	using _Alloc_traits = std::allocator_traits<_Alloc>;
	typename _Alloc_traits::template rebind_alloc<std::complex<_T>>
		_Cmpl_T_alloc;
	typename _Alloc_traits::template rebind_alloc<_T>
		_T_alloc;
	typename _Alloc_traits::template rebind_alloc<int>
		_int_alloc;
	*/
	using size_type = std::make_unsigned_t<int>;

	char const JOBZ   = (Z == nullptr) ? 'N' : 'V';
	char const RANGE  = 'A';
	char const UPLO   = 'U';
	_T   const VL     { 0.0 };
	_T   const VU     { 0.0 };
	int  const IL     = 0;
	int  const IU     = 0;
	_T   const ABSTOL { 0.0 }; //std::numeric_limits<T>::min();
	int        M      = 0;
	auto       ISUPPZ = 
		tcm::utils::_Storage<int, _Alloc>{2 * static_cast<size_type>(N)};
	// tcm::utils::allocate_workspace(_int_alloc, 2 * N);
	                 // std::make_unique<int[]>(2 * N);
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
			, ISUPPZ.data()
			, &_work_dummy, &LWORK
			, &_rwork_dummy, &LRWORK
			, &_iwork_dummy, &LIWORK
			, &INFO
		    );

		LWORK  = static_cast<int>(std::real(_work_dummy));
		LRWORK = static_cast<int>(_rwork_dummy);
		LIWORK = _iwork_dummy;
	}

	auto       WORK   = 
		tcm::utils::_Storage<std::complex<_T>, _Alloc>{
			static_cast<size_type>(LWORK) };
	// tcm::utils::allocate_workspace(_Cmpl_T_alloc, LWORK);
	                 // std::make_unique<std::complex<_T>[]>(LWORK);
	auto       RWORK  = 
		tcm::utils::_Storage<_T, _Alloc>{static_cast<size_type>(LRWORK)};
	// tcm::utils::allocate_workspace(_T_alloc, LRWORK);
	                 // std::make_unique<_T[]>(LRWORK);
	auto       IWORK  = 
		tcm::utils::_Storage<int, _Alloc>{static_cast<size_type>(LIWORK)};
	// tcm::utils::allocate_workspace(_int_alloc, LIWORK);
	                 // std::make_unique<int[]>(LIWORK);

	tcm::import::heevr<std::complex<_T>>
		( &JOBZ, &RANGE, &UPLO, &N
		, A, &LDA
		, &VL, &VU
		, &IL, &IU
		, &ABSTOL, &M
		, W, Z, &LDZ
		, ISUPPZ.data()
		, WORK.data(), &LWORK
		, RWORK.data(), &LRWORK
		, IWORK.data(), &LIWORK
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



template< class _T
        , class _Alloc = std::allocator<_T>
        >
auto heevr( std::size_t const n
          , _T* A, std::size_t const lda
          , utils::Base<_T>* W
          , _T* Z, std::size_t const ldz ) -> void
{
	MEASURE;

	heevr_impl<_Alloc>
		( boost::numeric_cast<int>(n)
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


/*
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
*/




#endif // TCM_HERMITIAN_HPP
