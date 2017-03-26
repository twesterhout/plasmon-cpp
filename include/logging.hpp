#ifndef TCM_LOGGING_HPP
#define TCM_LOGGING_HPP

#include <boost/date_time/posix_time/posix_time_types.hpp>

#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/smart_ptr/make_shared_object.hpp>

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/sync_frontend.hpp>
#include <boost/log/sinks/text_ostream_backend.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>

#include <boost/log/support/date_time.hpp>

namespace logging = boost::log;


namespace tcm {


using severity_level = boost::log::trivial::severity_level;


auto setup_console_logging()
{
	using namespace boost::log;
	using namespace boost::log::expressions;

	auto sink = add_console_log(std::clog);

    sink->set_formatter
	( stream 
		<< attr<unsigned int>("LineID")
		<< ": [" << format_date_time<boost::posix_time::ptime>
		            ("TimeStamp", "%Y-%m-%d %H:%M:%S")
		<< "] [" << trivial::severity
		<< "] " << smessage
	);

	add_common_attributes();
}



} // namespace tcm

#define LOG(lg, severity) BOOST_LOG_SEV(lg, tcm::severity_level::severity)



#endif // TCM_LOGGING_HPP
