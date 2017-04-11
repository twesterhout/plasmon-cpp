#ifndef TCM_BENCHMARK_HPP
#define TCM_BENCHMARK_HPP

#include <iomanip>
#include <chrono>
#include <unordered_map>
#include <thread>
#include <mutex>

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include <detail/config.hpp>


///////////////////////////////////////////////////////////////////////////////
/// \file benchmark.hpp
/// \brief Defines utilities related to benchmarking.
///
/// \detail To get an idea how long a certain part of a simulation took, we
/// measure the execution time of some functions. This can be turned on/off
/// by defining/undefining the #CONFIG_DO_MEASURE flag in detail/config.hpp 
/// file. Obtained measurements are saved in a global static table 
/// `tcm::timing::_global_impl_stats`. Please, __do not use__ this directly.
/// There is also an global static #std::mutex 
/// `tcm::timing::_global_impl_stats_mutex` protecting `_global_impl_stats`.
/// Please, __do not__ temper with it either.
///
/// There are two functions that should be used to manipulate this table.
/// * #update() is used to record benchmarks;
/// * #report() is used to report the results.
///
/// We also implement a #Timer class. It starts the timer at construction and
/// stops it and saves the results at destruction. To simplify creation of
/// #Timer objects a #TCM_MEASURE macro is provided.
///
/// __Example usage__:
/// \include src/timing_example.cpp
/// __Possible output__:
/// \code{.unparsed}
/// [-------------------------------------------------------]
/// [ obscure_namespace::foobarfoo() |              7.5e-08 ]
/// [                          foo() |           0.00513702 ]
/// [-------------------------------------------------------]
/// \endcode
///////////////////////////////////////////////////////////////////////////////


namespace tcm {

namespace timing {


static std::unordered_map< std::string
                         , std::chrono::duration<double> > _global_impl_stats;
static std::mutex _global_impl_stats_mutex;


///////////////////////////////////////////////////////////////////////////////
/// \brief Updates the global stats table in a thread-safe way.

/// \param func_name Name of the function how it will appear in the table.
/// \param delta_t   Extra seconds spent in \p func_name since the last call
///                  to #update(). \p delta_t is added to the old time. If this
///                  is the first call to #update() old time is taken to be 0.
/// 
/// This function has the following behavior:
/// \snippet include/benchmark.hpp Updating global table
///////////////////////////////////////////////////////////////////////////////
auto update( std::string func_name
           , std::chrono::duration<double> const delta_t ) -> void
{
	//! [Updating global table]
	std::lock_guard<std::mutex> lock{ _global_impl_stats_mutex };
	auto i = _global_impl_stats.emplace(std::move(func_name), 0.0).first;
	i->second += delta_t;
	//! [Updating global table]
}


///////////////////////////////////////////////////////////////////////////////
/// \brief Pretty-prints the global stats table to \p out.

/// \param out Output stream where to write to.
/// \warning #std::setw is used in the implementation which places some
///          constraints on _Stream.
///////////////////////////////////////////////////////////////////////////////
template <class _Stream>
auto report(_Stream& out) -> void
{
	std::lock_guard<std::mutex> lock{ _global_impl_stats_mutex };

	auto const get_length = [](auto const& s) noexcept { return s.size(); };
	auto const max_name_width = *boost::max_element( _global_impl_stats
		| boost::adaptors::map_keys
		| boost::adaptors::transformed(std::cref(get_length)) );
	auto constexpr max_time_width = 20;
	auto const hline = std::string(max_name_width + max_time_width + 5, '-');

	out << "[" << hline << "]\n";
	for (const auto& _element : _global_impl_stats) {
		out << "[ " << std::setw(max_name_width) << _element.first 
		    << " | " 
		    << std::setw(max_time_width) << _element.second.count() 
			<< " ]\n";
	}
	out << "[" << hline << "]\n";
}


struct Timer {
	Timer(std::string name)
	    : _name{ std::move(name) }
		, _start{ std::chrono::high_resolution_clock::now() }
		, _stop{}
	{}

	Timer(Timer const&) = delete;
	Timer& operator= (Timer const&) = delete;

	~Timer()
	{
		_stop = std::chrono::high_resolution_clock::now();
		update(std::move(_name), _stop - _start);
	}

private:
	using time_point = std::chrono::high_resolution_clock::time_point;
	
	std::string  _name;
	time_point   _start;
	time_point   _stop;
};


///////////////////////////////////////////////////////////////////////////////
/// \brief Simplifies the creation of #Timer objects.

/// If #CONFIG_DO_MEASURE is defined, creates a #Timer object, otherwise does
/// nothing. This allows to turn benchmarking on/off without changing source
/// files.
///////////////////////////////////////////////////////////////////////////////
#ifdef CONFIG_DO_MEASURE
#	define TCM_MEASURE(function_name) \
		tcm::timing::Timer _timer_temp_object_{function_name}
#else
#	define TCM_MEASURE(function_name) \
		do {} while(false)
#endif



} // namespace timing

} // namespace tcm


#endif // TCM_BENCHMARK_HPP
