#include <iostream>
#include <iterator>
#include <sstream>
#include <iomanip>
#include <map>

#include <boost/program_options.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/mpi.hpp>

#include <benchmark.hpp>
#include <logging.hpp>

#include <constants.hpp>
#include <lapack.hpp>
#include <matrix_serialization.hpp>
#include <dielectric_function_v2.hpp>


namespace po  = boost::program_options;
namespace mpi = boost::mpi;


///////////////////////////////////////////////////////////////////////////////
/// \brief Real field.
///////////////////////////////////////////////////////////////////////////////
using R = double;

///////////////////////////////////////////////////////////////////////////////
/// \brief Complex field.
///////////////////////////////////////////////////////////////////////////////
using C = std::complex<R>;


///////////////////////////////////////////////////////////////////////////////
/// \brief Returns the rank of the admin process.
///////////////////////////////////////////////////////////////////////////////
constexpr auto admin_rank() -> int { return 0; }



// ============================================================================
//                                   SETTINGS                                  
// ============================================================================

///////////////////////////////////////////////////////////////////////////////
/// \brief Class representing the configurations options.

/// It is a small wrapper around boost::program_options, which read the
/// options from command line during construction.
///////////////////////////////////////////////////////////////////////////////


auto init_options() -> po::options_description
{
	po::options_description description;
	description.add_options()
		( "help", "Produce the help message." )
		( "out.file.log" 
		, po::value<std::string>()->default_value("sample")
		, "Name of the file that will be used for logging. This is "
		  "not the full name, but rather a 'base'. The actual file name "
		  "will be \"[log-file].[PROCESS_RANK].log\"." )
		( "out.file.eps" 
		, po::value<std::string>()->default_value("Epsilon")
		, "Name of the file that will be used for saving the computed "
		  "dielectric functions. This is, again, a base rather than the "
		  "actual name. The actual file name will be "
		  "\"[eps-file].[PROCESS_RANK].bin\", which indicates that "
		  "dielectric functions will be stored in binary format of the "
		  "boost::serialization library." )
		( "in.file.energies"
		, po::value<std::string>()->required()
		, "Name of the BIN file where the eigenenergies of the hamiltonian"
		  " are read from. This file must be in the format of the "
		  "boost::serialization library. This option is REQUIRED." )
		( "in.file.states"
		, po::value<std::string>()->required()
		, "Name of the BIN file where the eigenstates of the hamiltonian "
		  "are read from. This file must be in the format of the "
		  "boost::serialization library. This option is REQUIRED." )
		( "in.file.potential"
		, po::value<std::string>()->required()
		, "Name of the BIN file where the interaction potential "
		  "is read from. This file must be in the format of the "
		  "boost::serialization library. This option is REQUIRED." )
		( "in.frequency.start"
		, po::value<R>()->required()
		, "Starting frequency in eV. Must be a real value. This option "
		  "is REQUIRED." )
		( "in.frequency.stop"
		, po::value<R>()->required()
		, "Stopping frequency in eV. Must be a real value. This option "
		  "is REQUIRED." )
		( "in.frequency.step"
		, po::value<R>()->required()
		, "Step in frequency in eV. Must be a real value." );
	description.add(tcm::init_constants_options<double>());
	return description;
}




template <class _Help, class _Proceed>
auto parse_command_line( int argc, char** argv
                       , _Help&& help
                       , _Proceed&& proceed ) -> bool
{
	auto const desc = init_options();
	po::variables_map vm;

	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

	if (vm.count("help")) {
		help(desc);
		return false;
	}

	po::notify(vm);
	proceed(vm);
	return true;
}




auto log_file_name(int const rank, std::string file_name_base)
{
	return file_name_base + "." + std::to_string(rank) + ".log";	
}

auto initialize_logging( int const rank
                       , std::string const& file_name_base) -> void
{
	using namespace boost::log;
	register_simple_formatter_factory<tcm::severity_level, char>("Severity");
	add_file_log
	( 
		keywords::file_name = log_file_name(rank, file_name_base),
		keywords::format = "%LineID%: [%TimeStamp%] [%Severity%] %Message%",
		keywords::auto_flush = true
	);

	add_common_attributes();
}


template <class _R, class _C, class _F>
struct IPackage {
	std::tuple<_R, _R, _R>                frequency_range;
	std::string                           log_file_name_base;
	std::string                           eps_file_name_base;
	tcm::Matrix<_F>                       E;
	tcm::Matrix<_C>                       Psi;
	tcm::Matrix<_C>                       V;
	std::map<std::string, _R>             constants;

private:
	friend boost::serialization::access;

	template<class _Archive>
	auto save(_Archive & ar, unsigned int const version) const -> void
	{
		ar << std::get<0>(frequency_range)
		   << std::get<1>(frequency_range)
		   << std::get<2>(frequency_range)
		   << log_file_name_base
		   << eps_file_name_base
		   << E 
		   << Psi
		   << V
		   << constants;
	}

	template<class _Archive>
	auto load(_Archive & ar, unsigned int const version) -> void
	{
		ar >> std::get<0>(frequency_range)
		   >> std::get<1>(frequency_range)
		   >> std::get<2>(frequency_range)
		   >> log_file_name_base
		   >> eps_file_name_base
		   >> E 
		   >> Psi
		   >> V
		   >> constants;
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()
};


template<class _T>
auto load_matrix(std::string const& file_name) -> tcm::Matrix<_T>
{
	std::ifstream in_stream{file_name};
	if(not in_stream)
		throw std::runtime_error{"Failed to open `" + file_name + "`."};
	boost::archive::binary_iarchive in_archive{in_stream};

	tcm::Matrix<_T> A;
	in_archive >> A;
	return A;
}




auto load_ipackage(po::variables_map const& vm) -> IPackage<R, C, R>
{
	return { std::make_tuple( vm["in.frequency.start"].as<R>()
	                        , vm["in.frequency.stop"].as<R>()
	                        , vm["in.frequency.step"].as<R>() )
	       , vm["out.file.log"].as<std::string>()
	       , vm["out.file.eps"].as<std::string>()
	       , load_matrix<R>(vm["in.file.energies"].as<std::string>())
	       , load_matrix<C>(vm["in.file.states"].as<std::string>())
	       , load_matrix<C>(vm["in.file.potential"].as<std::string>())
		   , tcm::load_constants<R, double, std::map<std::string, R>>(vm)
		   };
}


template <class _T, class _Logger>
auto cache( std::string const& message
          , tcm::Matrix<_T> const& X
		  , std::string const& file_name
          , _Logger & lg ) -> void
{
	LOG(lg, info) << "Caching: " << message << "...";

	std::ofstream out_stream{file_name};
	if(not out_stream)
		throw std::runtime_error{"Could not open `" + file_name + "`."};
	boost::archive::binary_oarchive out_archive{out_stream};
	out_archive << X;

	LOG(lg, info) << "Caching successfully finished.";
}


template<class _Real, class _Logger>
auto get_job( mpi::communicator const& world
            , std::tuple<_Real, _Real, _Real> const& range
			, _Logger & lg ) -> std::vector<_Real>
{	
	LOG(lg, info) << "Calculating homework...";

	auto const  rank  = world.rank();
	auto const  size  = world.size();
	auto const& begin = std::get<0>(range);
	auto const& end   = std::get<1>(range);
	auto const& step  = std::get<2>(range);

	std::vector<_Real> homework;
	for(auto i = 0; begin + i * step <= end; ++i) {
		if(i % size == rank) 
			homework.push_back(begin + i * step);
	}
	
	auto record = lg.open_record(boost::log::keywords::severity = 
	                                 tcm::severity_level::info);
	if (record) {
		boost::log::record_ostream stream{record};
		stream << "Need to perform calculations for the following "
			   << "frequencies: {";
		if(not homework.empty()) {
			for(std::size_t i = 0; i < homework.size() - 1; ++i)
				stream << homework[i] << ", ";
			stream << homework.back();
		}
		stream << "}";
		stream.flush();
		lg.push_record(std::move(record));
	}
	
	return homework;
}


template<class _Number, class _F, class _R, class _C, class _Logger>
auto calculate_single( _Number const omega
                     , tcm::Matrix<_F> const& E
					 , tcm::Matrix<_C> const& Psi
					 , tcm::Matrix<_C> const& V
					 , std::map<std::string, _R> const& cs
                     , _Logger & lg 
					 , std::string const& file_name_base ) -> void
{
	using namespace std::complex_literals;
	LOG(lg, info) << "Calculating dielectric function for omega = "
	              << omega << "...";

	auto const file_name_matrix = 
		file_name_base + "." + std::to_string(std::real(omega)) 
		+ ".matrix.bin";
	auto const file_name_eigenvalues = 
		file_name_base + "." + std::to_string(std::real(omega)) 
		+ ".eigenvalues.bin";
	auto const file_name_eigenstates = 
		file_name_base + "." + std::to_string(std::real(omega)) 
		+ ".eigenstates.bin";

	auto epsilon = tcm::dielectric_function::make(omega, E, Psi, V, cs, lg);
	cache("Dielectric function matrix", epsilon, file_name_matrix, lg);

	LOG(lg, info) << "Diagonalizing dielectric function for omega = "
	              << omega << "...";

	tcm::Matrix<_C> W{epsilon.height(), 1};
	tcm::Matrix<_C> Z{epsilon.height(), epsilon.height()};
	tcm::lapack::geev(epsilon, W, Z);

	cache("Dielectric function eigenvalues", W, file_name_eigenvalues, lg);
	cache("Dielectric function eigenstates", Z, file_name_eigenstates, lg);

	LOG(lg, info) << "Done for omega = " << omega << "!";
}



auto run( mpi::communicator & world
        , IPackage<R, C, R> & input ) -> void
{
	mpi::broadcast(world, input, admin_rank());

	initialize_logging(world.rank(), input.log_file_name_base);
	boost::log::sources::severity_logger<tcm::severity_level> lg;

	auto const homework = get_job<R>(world, input.frequency_range, lg);
	for(auto const& w : homework) {
		calculate_single( C{w, input.constants.at("tau")}
		                , input.E
						, input.Psi
						, input.V
						, input.constants
						, lg
						, input.eps_file_name_base );
	}

	auto record = lg.open_record(boost::log::keywords::severity = 
	                                 tcm::severity_level::info);
	if (record) {
		boost::log::record_ostream stream{record};
		stream << "Timings:\n";
		tcm::timing::report(stream);
		stream.flush();
		lg.push_record(std::move(record));
	}
}






int main(int argc, char** argv)
{
	mpi::environment env;
	mpi::communicator world;

	IPackage<R, C, R> input;
	bool proceed_with_calculation;
	if (world.rank() == admin_rank()) {
		proceed_with_calculation = 
			parse_command_line( argc, argv
			                  , [&env] (auto _desc) { std::cout << _desc << '\n'; }
			                  , [&input] (auto _vm) { input = load_ipackage(_vm); } );
	}
	
	mpi::broadcast(world, proceed_with_calculation, admin_rank());
	if(not proceed_with_calculation)
		return EXIT_SUCCESS;
	
	run(world, input);

	return EXIT_SUCCESS;
}
