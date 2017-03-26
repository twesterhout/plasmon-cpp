#ifndef TCM_LAPACK_WRAPPER_HPP
#define TCM_LAPACK_WRAPPER_HPP

#include <config.hpp>

#ifdef USING_INTEL_MKL
#	include <lapack_wrapper_mkl.hpp>
#else
#	ifdef USING_ATLAS
#		include <lapack_wrapper_atlas.hpp>
#	else
#		error "Need ATLAS or MKL"
#	endif
#endif


#endif // TCM_LAPACK_WRAPPER_HPP
