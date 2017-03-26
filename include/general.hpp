#ifndef TCM_GENERAL_HPP
#define TCM_GENERAL_HPP

#include <cassert>
#include <memory>
#include <algorithm>


#include <utils.hpp>
#include <lapack_wrapper.hpp>
#include <benchmark.hpp>





namespace tcm {


namespace lapack {


namespace {


template< class _T
        , class = std::enable_if_t
	      <   std::is_same<_T, float>() 
	       or std::is_same<_T, double>()
          >
        >
auto geev_impl( int const N
              , _T* A, int const LDA
	          , std::complex<_T>* W
			  , std::complex<_T>* VL, int const LDVL
			  , std::complex<_T>* VR, int const LDVR
	          , utils::Type2Type<_T> ) -> void
{
	if (N == 0) return;

	assert(N > 0);
	assert(A != nullptr and LDA >= N);
	assert(W != nullptr);
	assert(LDVL >= (VL == nullptr ? 1 : N));
	assert(LDVR >= (VR == nullptr ? 1 : N));

	auto       WR    = std::make_unique<_T[]>(N);
	auto       WI    = std::make_unique<_T[]>(N);
	char const JOBVL = VL == nullptr ? 'N' : 'V';
	char const JOBVR = VR == nullptr ? 'N' : 'V';
	int        LWORK = -1;
	int        INFO  = 0;

	{
		_T _work_dummy;

		tcm::import::geev<_T>
			( &JOBVL, &JOBVR, &N
			, A, &LDA
			, WR.get(), WI.get()
			, VL, &LDVL
			, VR, &LDVR
			, &_work_dummy, &LWORK
			, &INFO
			);

		LWORK = static_cast<int>(_work_dummy);
	}

	auto       WORK  = std::make_unique<_T[]>(LWORK);

	tcm::import::geev<_T>
	    ( &JOBVL, &JOBVR, &N
		, A, &LDA
		, WR.get(), WI.get()
		, VL, &LDVL
		, VR, &LDVR
		, WORK.get(), &LWORK
		, &INFO
		);

	if (INFO < 0) {
		throw std::invalid_argument{ "Argument #" + std::to_string(-INFO) 
		                           + " had an illegal value." };
	} 
	else if (INFO > 0) {
		throw std::runtime_error{"Call to ?GEEV failed."};
	}

	std::transform( WR.get(), WR.get() + N
	              , WI.get()
	              , W
	              , [](auto x, auto y){ return std::complex<_T>(x, y); }
	              );
}


template< class _T
        , class = std::enable_if_t
                  <   std::is_same<_T, float>() 
                   or std::is_same<_T, double>()
                  >
        >
auto geev_impl( int const N
              , std::complex<_T>* A, int const LDA
              , std::complex<_T>* W
              , std::complex<_T>* VL, int const LDVL
              , std::complex<_T>* VR, int const LDVR
	          , utils::Type2Type<std::complex<_T>> ) -> void
{
	if (N == 0) return;

	assert(N > 0);
	assert(A != nullptr and LDA >= N);
	assert(W != nullptr);
	assert(LDVL >= (VL == nullptr ? 1 : N));
	assert(LDVR >= (VR == nullptr ? 1 : N));

	char const JOBVL = VL == nullptr ? 'N' : 'V';
	char const JOBVR = VR == nullptr ? 'N' : 'V';
	auto       RWORK = std::make_unique<_T[]>(2 * N);
	int        LWORK = -1;
	int        INFO  = 0;

	{
		std::complex<_T> _work_dummy;

		tcm::import::geev<std::complex<_T>>
			( &JOBVL, &JOBVR, &N
			, A, &LDA
			, W
			, VL, &LDVL
			, VR, &LDVR
			, &_work_dummy, &LWORK
			, RWORK.get()
			, &INFO
		    );

		LWORK = static_cast<int>(std::real(_work_dummy));
	}

	auto       WORK  = std::make_unique<std::complex<_T>[]>(LWORK);

	tcm::import::geev<std::complex<_T>>
		( &JOBVL, &JOBVR, &N
		, A, &LDA
		, W
		, VL, &LDVL
		, VR, &LDVR
		, WORK.get(), &LWORK
		, RWORK.get()
	    , &INFO
	    );

	if (INFO < 0) {
		throw std::invalid_argument{ "Argument #" + std::to_string(-INFO) 
		                           + " had an illegal value." };
	} 
	else if (INFO > 0) {
		throw std::runtime_error{"Call to ?GEEV failed."};
	}
}

} // unnamed namespace


template<class _T>
inline
auto geev( std::size_t const n
         , _T* A, std::ptrdiff_t const lda
         , std::complex<utils::Base<_T>>* W
         , std::complex<utils::Base<_T>>* VL, std::ptrdiff_t const ldvl
         , std::complex<utils::Base<_T>>* VR, std::ptrdiff_t const ldvr
         ) -> void
{
	MEASURE;

	geev_impl( boost::numeric_cast<int>(n)
	         , A, boost::numeric_cast<int>(lda)
			 , W
			 , VL, boost::numeric_cast<int>(ldvl)
			 , VR, boost::numeric_cast<int>(ldvr)
			 , utils::Type2Type<_T>{} );
}



} // namespace lapack

} // namespace tcm



#endif // TCM_GENERAL_HPP
