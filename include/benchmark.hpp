#ifndef TCM_BENCHMARK_HPP
#define TCM_BENCHMARK_HPP

#include <iostream>
#include <chrono>
#include <unordered_map>
#include <thread>
#include <mutex>


///////////////////////////////////////////////////////////////////////////////
/// \file benchmark.hpp
/// \brief Defines utilities related to benchmarking.
///////////////////////////////////////////////////////////////////////////////


namespace tcm {

namespace timing {


///////////////////////////////////////////////////////////////////////////////
/// \brief Global object that represents statistics,
///        i.e. how long an execution of a particular function took.
///////////////////////////////////////////////////////////////////////////////
std::unordered_map< std::string, std::chrono::duration<double> > stats;
std::mutex stats_mutex;


///////////////////////////////////////////////////////////////////////////////
/// \brief Updates stats in a thread-safe way
///////////////////////////////////////////////////////////////////////////////
auto update( std::string const& func_name
           , std::chrono::duration<double> const delta_t ) -> void
{
	std::lock_guard<std::mutex> lock{stats_mutex};

	if (! stats.count(func_name) ) {
		stats[func_name] = std::chrono::duration<double>{ 0.0 };
	}

	stats[func_name] += delta_t;
}

///////////////////////////////////////////////////////////////////////////////
/// \brief Pretty-print stats to \p out.
///////////////////////////////////////////////////////////////////////////////
template<class _Stream>
auto report(_Stream& out) -> void
{
	std::lock_guard<std::mutex> lock{stats_mutex};

	out << "[******************************]\n"
	    << "[*] Benchmarking results: \n\n";
	for (const auto& _element : stats) {
		out << "[*] `" << _element.first << "` ----> " 
		    << _element.second.count() << " seconds.\n";
	}

	out << "[******************************]\n";
}


struct Measure {
	Measure(char const* name)
	    : _name{ name }
		, _start{ std::chrono::high_resolution_clock::now() }
		, _stop{}
	{}

	Measure(Measure const&) = delete;
	Measure& operator= (Measure const&) = delete;

	~Measure()
	{
		_stop = std::chrono::high_resolution_clock::now();
		update(_name, _stop - _start);
	}

private:
	char const*                                    _name;
	std::chrono::high_resolution_clock::time_point _start;
	std::chrono::high_resolution_clock::time_point _stop;
};


#ifdef DO_MEASURE
#	define MEASURE \
		tcm::timing::Measure _measure_it_temp_object_(__PRETTY_FUNCTION__)
#else
#	define MEASURE \
		do {} while(false)
#endif


} // namespace timing


} // namespace tcm


#endif // TCM_BENCHMARK_HPP
